#! /bin/bash

### This script aggregates the signal from a bigwig over a given binsize
### For that, either pass the binsize, or the bed file containing the bins

set -e # exit when any command fails
if [ $# -lt 3 ]; then echo 'Usage: aggregate_bigwig_over_windows.sh <in.bw> <chr.sizes> <windowSize> <optional: windows.bed>'; exit 1; fi

infile=$1
sizes=$2
windowSize=$3

if [ $# == 3 ]; then
	windowbed=$(basename "${sizes%.*}").${windowSize}bp.windows.bed
	if [ ! -f $windowbed ]; then
		# create bed file containing the bins
		bedtools makewindows -g $sizes -w $windowSize | awk -vOFS="\t" '{print $1,$2,$3,NR}' > $windowbed
	fi
else
	windowbed=$4
fi

x=$(basename "${infile%.*}")
outfile=${x}_averagedOver${windowSize}bp.bw

if [ -f $outfile ]; then
	echo $outfile exists
else
	echo Processing $x
	bigWigAverageOverBed $infile $windowbed ${x}_averagedOver${windowSize}bp.tab -bedOut=${x}_averagedOver${windowSize}bp.bed
	awk -F"\t" -vOFS="\t" '{print $1,$2,$3,$5}' ${x}_averagedOver${windowSize}bp.bed | sort -k1,1 -k2,2n > ${x}_averagedOver${windowSize}bp.sorted.bed
	bedGraphToBigWig ${x}_averagedOver${windowSize}bp.sorted.bed $sizes $outfile
	rm ${x}_averagedOver${windowSize}bp*.bed ${x}_averagedOver${windowSize}bp.tab
fi
