#! /bin/bash

### This script creates a .2bit and a .genome file either from a genome fasta file (if provided), or it tries to download it.

assembly_dir=$1
assembly_name=$2

if [ "$#" -eq 3 ]; then
	# from fasta
	echo compute 2bit from fasta
	fasta=$3
	faToTwoBit $fasta $assembly_dir/$assembly_name.2bit

elif [ "$#" -eq 2 ]; then
	# download
	url_UCSC=http://hgdownload.soe.ucsc.edu/goldenPath/$assembly_name/bigZips/$assembly_name.2bit
	url_UCSC_2=http://hgdownload.soe.ucsc.edu/gbdb/$assembly_name/$assembly_name.2bit
	url_damir=http://data.genereg.net/damir/$assembly_name/$assembly_name.2bit
	url_NHGRI=https://goldfish.nhgri.nih.gov/$assembly_name/$assembly_name/$assembly_name.2bit

	if [ ! -f $assembly_name.2bit ]; then
		echo "Download $assembly_name.2bit assembly"
		if wget -q --method=HEAD $url_UCSC; then
			url=$url_UCSC
		elif wget -q --method=HEAD $url_UCSC_2; then
			url=$url_UCSC_2
		elif wget -q --method=HEAD $url_damir; then
			url=$url_damir
		elif wget -q --method=HEAD $url_NHGRI; then
			url=$url_NHGRI
		fi
		if [[ $url != "" ]]; then
			wget -nc $url -O $assembly_dir/$assembly_name.2bit
		else
			echo ".2bit file not found on UCSC, NHGRI or Damir's repository."
			exit 1
		fi
	fi

else
	echo "Usage: ./get_assembly.sh <output_dir> <$assembly_name> <optional_genome_fasta_file>"
	exit 1
fi

# compute genome size
if [ -f $assembly_name.2bit ] && [ ! -f $assembly_name.sizes ]; then
	echo "Retrieving $assembly_name genome size"
	get_chromosome_sizes.R $assembly_dir/$assembly_name.2bit # script that creates tab-delimited .genome file with chromosomes and their lengths
else :
fi
