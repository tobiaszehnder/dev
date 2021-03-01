#! /bin/bash

[[ $# != 1 ]] && echo "Usage: ./get_chromosome_sizes.R <assembly.2bit>" && exit 1

twobit_file=$1
outfile=$(sed -e 's/\.2bit/\.sizes/g' <<< "$twobit_file")
twoBitInfo $twobit_file stdout | sort -k2rn > $outfile
