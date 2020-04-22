#! /bin/bash

### This script projects chunks of a GRB from zebrafish to mouse using Dijkstra's Shortest Path algorithm.
### Runtime: The projection of a ~ 1 MB GRB binned into 5kb bins running on 30 threads takes about 25 minutes.

##### TO DO: implement choice of ref / qry species

[[ $# != 4 ]] && { echo "Usage: ./project_grb.sh grb_name grb_coords.bed binsize nthreads"; exit 1; }

# parse args
grb_name=$1
grb_coords=( $(head -1 $2) ) # old approach (direct coordiniates): IFS='[:-]'; grb_coords=($2); unset IFS; # separate chrN:xxx-yyy and store in array ( chrN xxx yyy )
binsize=$3
nthreads=$4

# bin GRB
grb_chunks=(`seq ${grb_coords[1]} $binsize ${grb_coords[2]}`)
l="${#grb_chunks[@]}"
[[ $nthreads -gt $l ]] && nthreads=$l

# project GRB bins (parallelize)
mkdir -p tmp
ERT=$((4*l/nthreads)) # based on a estimated runtime of 4 minutes per job
echo "Projecting $2 in bins of size $3 using $4 threads in parallel"
echo "Estimated runtime: $((ERT/60))h $(bc <<< $ERT%60)min"
for j in `seq 0 $nthreads $l`; do
	for i in `seq $j $((j+nthreads-1))`; do
		[[ $i -lt $l ]] || break
		id="${grb_name}_${i}"
		coord=${grb_coords[0]}:${grb_chunks[$i]}
		echo $id $coord
		project_dijkstra.py $coord $id &
	done
	wait
done

# concatenate tmp output files to one file, delete tmp files
outfile=${grb_name}.proj
head -2 "tmp/${grb_name}_1.proj.tmp" > $outfile;  eval tail -n 1 -q "tmp/${grb_name}_{0..$((l-1))}.proj.tmp" >> $outfile
rm -r tmp
echo "Done"
