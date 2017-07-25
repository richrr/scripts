#!/bin/bash
#sendsketch in=<infile> out=<outfile>

usage(){
echo "
Written by Brian Bushnell
Last modified  February 27, 2017

Description:  Compares a Sketches to others, and prints their kmer identity.
The files can be sketches made by sketch.sh, or fasta files.
This is similar to comparesketch.sh, but connects to a remote server
which is holding the reference sketeches in memory.

Usage:  sendsketch.sh in=file


Standard parameters:
in=file             Sketch or fasta file to compare.
size=4000           Size of sketch to generate.
amino=f             Use amino acid mode.
mode=single         Possible modes, for fasta input:
                       single: Generate one sketch per file.
                       sequence: Generate one sketch per sequence.
address=            Address of remote server.  Default adddress:
                    https://taxonomy.jgi-psf.org/sketch

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
	freeRam 3200m 84
	z="-Xmx${RAM}m"
	z2="-Xms${RAM}m"
}
calcXmx "$@"

sendsketch() {
	if [[ $NERSC_HOST == genepool ]]; then
		module unload oracle-jdk
		module load oracle-jdk/1.8_64bit
	fi
	local CMD="java $EA $z -cp $CP sketch.SendSketch $@"
	echo $CMD >&2
	eval $CMD
}

sendsketch "$@"
