#!/bin/bash
#https://askubuntu.com/questions/693409/how-can-i-extract-multiple-gzip-files-in-directory-and-subdirectories

# usage:  bash /nfs3/PHARM/Morgun_Lab/richrr/scripts/bash/uncompress-files.sh /nfs3/PHARM/Morgun_Lab/richrr/Milena_data_050917/test_few_samples/test_samps/ /nfs3/PHARM/Morgun_Lab/richrr/Milena_data_050917/test_few_samples/yuka_script_results/
# make sure that the outfolder has "/" at the end

INFOLDER=$1 
OUTFOLDER=$2

cd $INFOLDER

if [[ "$OUTFOLDER" == */ ]]; then
    echo "outfolder name ok"
else
    OUTFOLDER="$OUTFOLDER/"
fi


if [ ! -d "$OUTFOLDER" ]; then
  # Control will enter here if $DIRECTORY doesn't exist.
  mkdir $OUTFOLDER
fi

#for f in *.gz; do gunzip -c $f > $OUTFOLDER"${f%.*}" ; done

for f in *.gz; do
 gunzip -c $f > $OUTFOLDER"${f%.*}"
 #echo "$OUTFOLDER$f"
done
