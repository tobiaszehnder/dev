#! /bin/bash

### This script works analog to propagate_anchors_wrapper.sh and project_grb.sh to use SAPP for a bed of multiple GRBs that all get binned to binsize and then propagated
[[ $# != 6 ]] && { echo "Usage: ./propagate_anchors_grbs.sh grb_coords.bed reference_species query_species binsize path_pwaln_pkl nthreads"; exit 1; }

# parse args
grb_bed=$1
ref=$2
qry=$3
binsize=$4
path_pwaln_pkl=$5

# check if bed-file contains at least 4 columns
[ $(awk '{print NF}' $grb_bed | sort -nu | head -n 1) -lt 4 ] && echo "Error: bed-file must contain at least 4 columns with the 4th being the GRB ID / name." && exit 1

# loop through bed-file
while IFS='' read -r LINE || [ -n "${LINE}" ]; do
	nthreads=$6 # reset to the original input if the variable was set to $l in the previous iteration
	# put coordinates of current line into bash array
	grb_coords=($LINE)
	grb_name=${grb_coords[3]}
	outfile=${grb_name}_$((binsize/1000))kb.aprop
	[ -f $outfile ] && echo "$outfile exists already, skip." && continue
	
	# bin GRB
	grb_chunks=(`seq $((grb_coords[1]/binsize*binsize)) $binsize $(((grb_coords[2]/binsize+1)*binsize+1))`)
	l="${#grb_chunks[@]}"
	[ $nthreads -gt $l ] && nthreads=$l

	# project GRB bins (parallelize)
	mkdir -p tmp
	ERT=$(printf "%.0f" "$(echo "30*$l/$nthreads" | bc -l)") # based on a estimated average runtime of 30 sec per job
	echo "Propagating anchors of $grb_name from $ref to $qry in $l bins of size $binsize using $nthreads threads in parallel"
	echo "Estimated runtime: $((ERT/3600))h $(bc <<< $((ERT/60)))m $(bc <<< $ERT%60)s"
	sem_id="project_${grb_name}_$(hostname)_${RANDOM}"
	starttime=$(date -u '+%s')
	for i in `seq 0 $((l-1))`; do
		id="${grb_name}_${i}"
		coord=${grb_coords[0]}:${grb_chunks[$i]}
		echo $id $coord
		sem --id $sem_id -j${nthreads} --timeout 80 propagate_anchors.py $ref $qry $coord $id $path_pwaln_pkl tmp # sem is an alias for parallel --semaphore. A counting semaphore will allow a given number of jobs to be started in the background.
	done
	sem --id $sem_id --wait # wait until all sem jobs are completed before moving on
	endtime=$(date -u '+%s')
	difftime=$(date -u --date @$((endtime-starttime)) '+%-Hh %-Mm %-Ss')
	echo "Effective runtime: ${difftime}"

	# concatenate tmp output files to one file, delete tmp files
	head -3 "tmp/${grb_name}_0.aprop" > $outfile
	for file in `ls tmp/${grb_name}_*.aprop | sort -V`; do eval tail -n+4 -q $file >> $outfile; done 
	# rm -r tmp
	echo "Done"
done < $grb_bed
