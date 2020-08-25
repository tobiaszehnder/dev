#! /bin/bash

# This script takes a bed file of point coordinates and projects them from a reference to a query.
# Make sure path_pwaln_pkl contains the desired species set.

[[ $# != 6 ]] && { echo "Usage: ./project_regions.sh coords.bed reference_species query_species half_life_distance path_pwaln_pkl nthreads"; exit 0; }

bed=$1
ref=$2
qry=$3
half_life_distance=$4
path_pwaln_pkl=$5
nthreads=$6

# check if bed-file contains at least 4 columns
[ $(awk '{print NF}' $bed | sort -nu | head -n 1) -lt 4 ] && echo "Error: bed-file must contain at least 4 columns with the 4th being the GRB ID / name." && exit 0

mkdir -p tmp
l=$(< $bed wc -l)
ERT=$(printf "%.0f" "$(echo "1*$l/$nthreads" | bc)") # based on a estimated average runtime of 1 min per job
echo "Estimated runtime: $((ERT/60))h $(bc <<< $ERT%60)m"
sem_id="project_coordinates_$(hostname)_${RANDOM}"
starttime=$(date -u '+%s')

# loop through bed-file
while IFS='' read -r LINE || [ -n "${LINE}" ]; do
	bed_row=($LINE)
	id=${bed_row[3]}
	coord=${bed_row[0]}:$(($((${bed_row[1]}+${bed_row[2]}))/2)) # center of region
	echo $id $coord
	sem --id $sem_id -j${nthreads} project_dijkstra.py $ref $qry $coord $id $half_life_distance $path_pwaln_pkl # sem is an alias for parallel --semaphore. A counting semaphore will allow a given number of jobs to be started in the background.
done < $bed

sem --id $sem_id --wait # wait until all sem jobs are completed before moving on
endtime=$(date -u '+%s')
difftime=$(date -u --date @$((endtime-starttime)) '+%-Hh %-Mm %-Ss')
echo "Effective runtime: ${difftime}"

# concatenate tmp output files to one file, delete tmp files
outfile=${bed/bed/proj}
ids=($(cut -f4 $bed))
head -3 "tmp/${ids[0]}.proj.tmp" > $outfile
for id in ${ids[@]:1}; do eval tail -n 1 -q "tmp/${id}.proj.tmp" >> $outfile; done
rm -r tmp
echo "Done"
