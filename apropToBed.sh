#! /bin/bash

### This script converts the genomic coordinates of the anchor points from the SAPP.aprop output to a BED format.
### Pass the .aprop file and the species of interest ("reference" or "target")

# parse arguments
[[ $# != 2 ]] && echo "Usage: apropToBed.sh infile.aprop species" && exit 1
infile=$1
species=$2
outfile=$(basename "${infile%.*}")_anchors_$species.bed

awk -F'\t' -vOFS='\t' -v x=$species 'BEGIN{c=1}{a=c;if(a==2){if($7==x){if($10=="+"){print $8,b-1,$9,$1,0,$10;c=1} else {print $8,$9-1,b,$1,0,$10;c=1}}} if(a==1){if($7==x){b=$9;c=2}}}' $infile > $outfile
echo "Extracted ${species} anchors and written to ${outfile}"
