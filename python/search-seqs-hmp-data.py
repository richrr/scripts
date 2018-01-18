import sys
import os
from utils import * 
import subprocess


# [Linux@samwise HMP]$ SGE_Batch -c "python ~/Morgun_Lab/richrr/scripts/python/search-seqs-hmp-data.py v3v5post-fornix-lmd.txt" -m 100G -F 200G -M rodrrich@oregonstate.edu -r sge_search_seqs_log

patterns = read_file(sys.argv[1])

reslist = list()
reslist.append("pattern\tseqs")

seqlist = list()

for p in patterns:
    p = p.strip()
    cmd = 'grep -c ' + p + ' seqs_v35.fna'
    #print cmd
    res = '0'
    # https://stackoverflow.com/questions/20983498/subprocess-check-output-with-grep-command-fails-when-grep-finds-no-matches
    try:
        res = subprocess.check_output(cmd, shell=True)
        res = res.strip()
    except subprocess.CalledProcessError as e:
        if e.returncode > 1:
            raise
    #print res
    out = p + '\t' + res
    print out
    reslist.append(out)


    cmd = 'grep -A1 ' + p + ' seqs_v35.fna'
    #print cmd
    res = '0'
    try:
        res = subprocess.check_output(cmd, shell=True)
        res = res.strip()
        seqlist.append(res)
    except subprocess.CalledProcessError as e:
        if e.returncode > 1:
            raise
    print res


writeLIST_to_file(reslist, "pattern-read-counts.tsv")
writeLIST_to_file(seqlist, "pattern-reads.txt")
