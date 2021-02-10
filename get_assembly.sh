#! /bin/bash

### This script creates a .2bit and a .genome file either from a genome fasta file (if provided), or it tries to download it.
### Usage: ./get_assembly.sh <output_dir> <assembly_name> <optional_genome_fasta_file>

assembly_dir=$1
assembly_name=$2

if [[ ! -e $assembly_dir/$assembly_name.2bit ]]; then
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
			if wget -q --method=HEAD $url_UCSC; then wget -nc $url_UCSC -O $assembly_dir/$assembly_name.2bit;
			elif wget -q --method=HEAD $url_UCSC_2; then wget -nc $url_UCSC_2 -O $assembly_dir/$assembly_name.2bit;
			elif wget -q --method=HEAD $url_damir; then wget -nc $url_damir -O $assembly_dir/$assembly_name.2bit;
			elif wget -q --method=HEAD $url_NHGRI; then wget -nc $url_NHGRI -O $assembly_dir/$assembly_name.2bit;
	   		else echo ".2bit file not found on UCSC, NHGRI or Damir's repository." && exit 1;
			fi
		fi
	else
		echo "Usage: ./get_assembly.sh <output_dir> <assembly_name> <optional_genome_fasta_file>"
		exit 1
	fi
else
	echo "File exists: $assembly_dir/$assembly_name.2bit"
fi

# compute genome size
if [ -f $assembly_dir/$assembly_name.2bit ]; then
	if [ ! -f $assembly_dir/$assembly_name.sizes ]; then
		echo "Retrieving $assembly_name genome size"
		get_chromosome_sizes.R $assembly_dir/$assembly_name.2bit # script that creates tab-delimited .genome file with chromosomes and their lengths
	else
		echo "File exists: $assembly_dir/$assembly_name.sizes"
	fi
fi
