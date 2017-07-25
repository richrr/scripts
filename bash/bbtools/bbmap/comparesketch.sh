#!/bin/bash
#comparesketch in=<infile> out=<outfile>

usage(){
echo "
Written by Brian Bushnell
Last modified  February 27, 2017

Description:  Compares a Sketches to others, and prints their kmer identity.
The files can be sketches made by sketch.sh, or fasta files.
It's recommended to first sketch references with sketch.sh for large files,
or when taxonomic information is desired.


Usage:  comparesketch.sh in=<file,file,file...> ref=<file,file,file...>
Alternative:  comparesketch.sh in=<file,file,file...> file file file
Alternative:  comparesketch.sh in=<file,file,file...> *.sketch.gz


Standard parameters:
in=<file,file...>   List of sketches to compare.
ref=<file,file...>  List of sketches to compare against.  Files given without
                    a prefix (ref=) will be treated as references,
                    so you can use *.sketch.gz without ref=.
threads=auto        Use this many threads for comparison.
index=auto          Index the sketches for much faster searching.
                    Requires more memory and adds startup time.
                    Recommended to be turned on for comparing many sketches.
amino=f             Use amino acid mode.
minhits=3           (hits) Only report records with at least this many hits.
minid=0.002         (id) Only report records with at least this identity (0-1).
records=100         Report at most this many best-matching records.
mode=single         Possible modes, for fasta input:
                       single: Generate one sketch per file.
                       sequence: Generate one sketch per sequence.

Taxonomy-related parameters:
tree=<file>         Specify a TaxTree file.  On Genepool, use tree=auto.
                    Only necessary for print results by taxonomy.
                    Assumes comparisons are done against reference sketches
                    with known taxonomy information.
printtaxa=t         Print full taxonomy of each record.
format=0                0: Print taxonomy as one line per level.
                        1: Print taxonomy as one line, semicolon-delimited.
level=2             Only report the best record per taxa at this level:
                       -1: disabled
                        1: subspecies
                        2: species
                        3: genus
                       ...etc

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

comparesketch() {
	if [[ $NERSC_HOST == genepool ]]; then
		module unload oracle-jdk
		module load oracle-jdk/1.8_64bit
	fi
	local CMD="java $EA $z -cp $CP sketch.CompareSketch $@"
	echo $CMD >&2
	eval $CMD
}

comparesketch "$@"
