#!/bin/bash
#plotgc in=<file> out=<file> ref=<ref file>

usage(){
echo "
Written by Brian Bushnell
Last modified February 27, 2017

Description:  Prints sequence gc content once per interval.

Usage:  plotgc.sh in=<input file> out=<output file>


Parameters:
in=<file>           Input file. in=stdin.fa will pipe from stdin.
out=<file>          Output file.  out=stdout will pipe to stdout.
interval=1000       Interval length.
offset=0            Position offset.  For 1-based indexing use offset=1.
printshortbins=t    (psb) Print gc content for the last bin of a contig
                    even when shorter than interval.

Java Parameters:

-Xmx                This will be passed to Java to set memory usage, overriding 
                    the program's automatic memory detection. -Xmx20g will 
                    specify 20 gigs of RAM, and -Xmx200m will specify 200 megs.  
                    The max is typically 85% of physical memory.

There is a changelog at /bbmap/docs/changelog_plotgc.txt
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
NATIVELIBDIR="$DIR""jni/"

z="-Xmx1g"
z2="-Xms1g"
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
	freeRam 1400m 42
	z="-Xmx${RAM}m"
	z2="-Xms${RAM}m"
}
calcXmx "$@"

plotgc() {
	if [[ $NERSC_HOST == genepool ]]; then
		module unload oracle-jdk
		module load oracle-jdk/1.8_64bit
		module load pigz
	fi
	local CMD="java $EA $z -cp $CP driver.PlotGC $@"
	echo $CMD >&2
	eval $CMD
}

plotgc "$@"
