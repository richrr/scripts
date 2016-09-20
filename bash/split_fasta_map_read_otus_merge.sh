mkdir uc_files

PARTS=$1
FILEN=$2

# split the fasta file into parts
perl ~/scripts/perl/fasta-splitter.pl -n-parts $PARTS $FILEN

# Map reads (including singletons) back to OTUs
for i in $(ls *.fasta)
   do 
   ~/bin/usearch61 -usearch_global $i -db ../otus97_filtered.fa -strand plus -id 0.97 -uc uc_files/$i.map.uc
   done

#http://stackoverflow.com/questions/13210880/replace-one-substring-for-another-string-in-shell-script
OLD='../'
NEW=''
OUTFILEN=${FILEN//$OLD/$NEW}
# merge the output files into one
cat uc_files/*.map.uc > ../$OUTFILEN.map.uc
