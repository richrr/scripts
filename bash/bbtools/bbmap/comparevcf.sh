#!/bin/bash
#comparevcf in=<file.sam> out=<file.vcf>

usage(){
echo "
Written by Brian Bushnell
Last modified January 20, 2017

Description:  Performs set operations on VCF files:
Union, intersection, and subtraction.

Usage:  comparevcf.sh in=<file,file,...> out=<file>


I/O parameters:
in=<file>       Input; must be at least 2 files.
out=<file>      Output file.
shist=<file>    (scorehist) Output for variant score histogram.
overwrite=f     (ow) Set to false to force the program to abort rather than
bgzip=f         Use bgzip for gzip compression.

Processing Parameters (choose one only):
subtract=t      Subtract all other files from the first file.
union=f         Make a union of all files.
intersection=f  Make an intersection of all files.
addsamples=t    Include all samples in the output lines. (TODO)


Java Parameters:
-Xmx            This will be passed to Java to set memory usage, overriding the program's automatic memory detection.
                -Xmx20g will specify 20 gigs of RAM, and -Xmx200m will specify 200 megs.  The max is typically 85% of physical memory.

Please contact Brian Bushnell at bbushnell@lbl.gov if you encounter any problems.
"
}

pushd . > /dev/null
DIR="${BASH_SOURCE[0]}"
while [ -h "$DIR" ]; do
  cd "$(dirname "$DIR")"
  DIR="$(readlink "$(basename "$DIR")")"
done
cd "$(dirname "$DIR")"
DIR="$(pwd)/"
popd > /dev/null

#DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )/"
CP="$DIR""current/"

z="-Xmx4g"
z2="-Xms4g"
EA="-ea"
set=0

if [ -z "$1" ] || [[ $1 == -h ]] || [[ $1 == --help ]]; then
	usage
	exit
fi

calcXmx () {
	source "$DIR""/calcmem.sh"
	parseXmx "$@"
	if [[ $set == 1 ]]; then
		return
	fi
	freeRam 4000m 84
	z="-Xmx${RAM}m"
	z2="-Xms${RAM}m"
}
calcXmx "$@"

comparevcf() {
	if [[ $NERSC_HOST == genepool ]]; then
		module unload oracle-jdk
		module load oracle-jdk/1.8_64bit
		module load pigz
	fi
	local CMD="java $EA $z $z2 -cp $CP var2.CompareVCF $@"
	echo $CMD >&2
	eval $CMD
}

comparevcf "$@"
