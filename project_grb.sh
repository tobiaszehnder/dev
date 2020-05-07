#! /bin/bash

### This script projects chunks of a GRB from zebrafish to mouse using Dijkstra's Shortest Path algorithm.
### Runtime: The projection of a ~ 1 MB GRB binned into 2kb bins running on 30 threads takes about 1 hour.
### New: you can pass multiple GRB coordinates in a bed-file. Make sure that the 4th column of the bed-file contains an ID / name.

[[ $# != 5 ]] && { echo "Usage: ./project_grb.sh grb_coords.bed reference_species query_species binsize nthreads"; exit 1; }

# parse args
grb_bed=$1
ref=$2
qry=$3
binsize=$4
nthreads=$5

# check if bed-file contains at least 4 columns
[ $(awk '{print NF}' $grb_bed | sort -nu | head -n 1) -lt 4 ] && echo "Error: bed-file must contain at least 4 columns with the 4th being the GRB ID / name." && exit 1

# loop through bed-file
while IFS='' read -r LINE || [ -n "${LINE}" ]; do
	# put coordinates of current line into bash array
	grb_coords=($LINE)
	grb_name=${grb_coords[3]}
	outfile=${grb_name}_$((binsize/1000))kb.proj
	[ -f $outfile ] && echo "$outfile exists already, skip." && continue
	
	# bin GRB
	grb_chunks=(`seq ${grb_coords[1]} $binsize ${grb_coords[2]}`)
	l="${#grb_chunks[@]}"
	[ $nthreads -gt $l ] && nthreads=$l

	# project GRB bins (parallelize)
	mkdir -p tmp
	ERT=$((3*l/nthreads)) # based on a estimated runtime of 3 minutes per job
	echo "Projecting $grb_name in $l bins of size $binsize using $nthreads threads in parallel"
	echo "Estimated runtime: $((ERT/60))h $(bc <<< $ERT%60)min"
	sem_id="project_${grb_name}" # $(hostname)_$(tty)" 
	for i in `seq 0 $((l-1))`; do
		id="${grb_name}_${i}"
		coord=${grb_coords[0]}:${grb_chunks[$i]}
		echo $id $coord
		sem --id $sem_id -j${nthreads} project_dijkstra.py $ref $qry $coord $id # sem is an alias for parallel --semaphore. A counting semaphore will allow a given number of jobs to be started in the background.
	done
	sem --id $sem_id --wait # wait until all sem jobs are completed before moving on

	# concatenate tmp output files to one file, delete tmp files
	head -2 "tmp/${grb_name}_0.proj.tmp" > $outfile;  eval tail -n 1 -q "tmp/${grb_name}_{0..$((l-1))}.proj.tmp" >> $outfile
	rm -r tmp
	echo "Done"
done < $grb_bed
