#! /bin/bash

# This script takes a bed file of point coordinates and propagates their anchors such that the anchor span in the target species is minimized.
# Make sure path_pwaln_pkl contains the desired species set.

[[ $# != 5 ]] && { echo "Usage: ./propagate_anchors_wrapper.sh coords.bed reference_species target_species path_pwaln_pkl nthreads"; exit 0; }

bed=$1
reference=$2
target=$3
path_pwaln_pkl=$4
nthreads=$5

# check if bed-file contains at least 4 columns
[ $(awk '{print NF}' $bed | sort -nu | head -n 1) -lt 4 ] && echo "Error: bed-file must contain at least 4 columns with the 4th being the GRB ID / name." && exit 0

mkdir -p tmp
l=$(< $bed wc -l)
ERT=$(printf "%.0f" "$(echo "30*$l/$nthreads" | bc -l)") # based on a estimated average runtime of 30 sec per job
echo "Estimated runtime: $((ERT/3600))h $(bc <<< $((ERT/60)))m $(bc <<< $ERT%60)s"
sem_id="propagate_anchors_$(hostname)_${RANDOM}"
starttime=$(date -u '+%s')

# loop through bed-file
while IFS='' read -r LINE || [ -n "${LINE}" ]; do
	bed_row=($LINE)
	id=${bed_row[3]}
	coord=${bed_row[0]}:$(($((${bed_row[1]}+${bed_row[2]}))/2)) # center of region
	echo $id $coord
	sem --id $sem_id -j${nthreads} --timeout 100 propagate_anchors.py $reference $target $coord $id $path_pwaln_pkl tmp
done < $bed

sem --id $sem_id --wait # wait until all sem jobs are completed before moving on
endtime=$(date -u '+%s')
difftime=$(date -u --date @$((endtime-starttime)) '+%-Hh %-Mm %-Ss')
echo "Effective runtime: ${difftime}"

# concatenate tmp output files to one file, delete tmp files
outfile=${bed/bed/aprop}
ids=($(cut -f4 $bed))
cat "tmp/${ids[0]}.aprop" > $outfile
for id in ${ids[@]:1}; do eval tail -n+4 -q "tmp/${id}.aprop" >> $outfile; done
rm -r tmp
echo "Done"
