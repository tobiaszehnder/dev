#! /bin/bash

if [ "$#" -ne 1 ]; then
    echo "Usage: ./get_assembly.sh <assembly_name>"
	exit 1
fi

assembly_dir=/project/wig/tobias/reg_evo/data/assembly
url_UCSC=http://hgdownload.soe.ucsc.edu/goldenPath/$1/bigZips/$1.2bit
url_UCSC_2=http://hgdownload.soe.ucsc.edu/gbdb/$1/$1.2bit
url_damir=http://data.genereg.net/damir/$1/$1.2bit
url_NHGRI=https://goldfish.nhgri.nih.gov/$1/$1/$1.2bit
echo $url_UCSC_2
# download .2bit file
if [ ! -f $1.2bit ]; then
	echo "Downloading $1.2bit assembly"
	if wget -q --method=HEAD $url_UCSC; then
		url=$url_UCSC
	elif wget -q --method=HEAD $url_UCSC_2; then
		url=$url_UCSC_2
	elif wget -q --method=HEAD $url_damir; then
		url=$url_damir
	elif wget -q --method=HEAD $url_NHGRI; then
		url=$url_NHGRI
	fi
	if [ $url != "" ]; then
		wget -nc $url -O $assembly_dir/$1.2bit
	else
		echo ".2bit file not found on UCSC, NHGRI or Damir's repository."
	fi
fi

# compute genome size
if [ -f $1.2bit ] && [ ! -f $1.sizes ]; then
	echo "Retrieving $1 genome size"
	get_chromosome_sizes.R $assembly_dir/$1.2bit # script that creates tab-delimited .genome file with chromosomes and their lengths
else :
fi
