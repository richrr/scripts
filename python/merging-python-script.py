import sys
import os
import pandas as pd
import multiprocessing as mp
import csv


# this code is written for the merged file with combined pval & fdr. although it could have been written for the file without comb fisher and fdr,
# it is easier to have the output with the comb pval and fdr and use what we need rather than have to search them in the merged file with comb pval and fdr
# or run the next (create network) command to calc the combined pval and fdr.


# cd /nfs3/PHARM/Morgun_Lab/richrr/Cervical_Cancer/analysis/merged/corr/gexpress/stage-ltest_corr/p1

"""
SGE_Batch -c "python ~/Morgun_Lab/richrr/scripts/python/merging-python-script.py merged_gexp_sp_corr_p1_FolChMedian_merged-parallel-output.csv-comb-pval-output.csv 1 2" -m 150G -F 100G -r log_merge-py_1 -q biomed -M rodrrich@oregonstate.edu -P 8

SGE_Batch -c "python ~/Morgun_Lab/richrr/scripts/python/merging-python-script.py merged_gexp_sp_corr_p1_FolChMedian_merged-parallel-output.csv-comb-pval-output.csv 2 2" -m 150G -F 100G -r log_merge-py_2 -q biomed -M rodrrich@oregonstate.edu -P 8

SGE_Batch -c "python ~/Morgun_Lab/richrr/scripts/python/merging-python-script.py merged_gexp_sp_corr_p1_FolChMedian_merged-parallel-output.csv-comb-pval-output.csv 3 2" -m 150G -F 100G -r log_merge-py_3 -q biomed -M rodrrich@oregonstate.edu -P 8

SGE_Batch -c "python ~/Morgun_Lab/richrr/scripts/python/merging-python-script.py merged_gexp_sp_corr_p1_FolChMedian_merged-parallel-output.csv-comb-pval-output.csv 4 2" -m 150G -F 100G -r log_merge-py_4 -q biomed -M rodrrich@oregonstate.edu -P 8

"""


infile = sys.argv[1]
analysis = "Analys " + sys.argv[2] + " "
numb_datasets = int(sys.argv[3])


# get the header line form the big file and decide which (analysis) columns to use
header_line = ''
with open(infile, 'r') as f:
    header_line = f.readline().strip()

selcted_cols = [i for i, s in enumerate(header_line.split(',')) if analysis in s]  #[s for s in header_line.split(',') if analysis in s]
# get the lowest and highest and make range out of it
# this way you get the combinedpval and combined fdr cols
selcted_cols = range(min(selcted_cols), max(selcted_cols)+1)
selcted_cols.insert(0, 0) # explicitly adding the row id cols

print selcted_cols

header_for_print = [header_line.split(',')[i] for i in selcted_cols]
print header_for_print


def process(df):
    res = list()
    for row in df.itertuples():
        #print row
        
        corrs = row[1:numb_datasets+1]
        corrs_flag = 0
                
        # write some condition to check for NA
        pos = sum(float(num) > 0 for num in corrs)
        neg = sum(float(num) < 0 for num in corrs)
        
        #print pos, neg
        if len(corrs) == pos and not len(corrs) == neg:
            #print "pos"
            corrs_flag = 1
        if len(corrs) == neg and not len(corrs) == pos:
            #print "neg"
            corrs_flag = 1

        if corrs_flag == 1: 
            res.append(row)
        
    return res

counter=0
pool = mp.Pool(30) # use 30 processes
funclist = []

# http://gouthamanbalaraman.com/blog/distributed-processing-pandas.html
#for chunck_df in pd.read_csv(infile, chunksize=100, usecols=range(5), index_col=0):
for chunck_df in pd.read_csv(infile, chunksize=100000, usecols=selcted_cols, index_col=0):
    counter = counter + 1
    print counter
    #print chunck_df
    # process each data frame
    f = pool.apply_async(process,[chunck_df])
    funclist.append(f)

#result = list()

OUTfile = infile + analysis.replace(" ", "_") + '-same-dir-corrs.csv'

with open(OUTfile, 'w') as of:
    writer = csv.writer(of, delimiter=',', lineterminator='\n')
    writer.writerow(header_for_print)
    for f in funclist:
        csvd = f.get(timeout=10000) # timeout in 10000 seconds
        #result.extend(csvd) 
        writer.writerows(csvd)


#print result

# quick and dirty command to get the first column of the file:
cutcmd = "cut -d, -f 1 " + OUTfile + " > " + OUTfile + "-ids.csv"
os.system(cutcmd)  


print "Done"


""" # sequential
corrs_dict = dict() # satisfies corr direction
counter = 0
# with open(in_filename) as in_f, open(out_filename, 'w') as out_f
with open(infile) as f: 
    for line in f:
        counter = counter + 1
        line = line.strip()
        print line
        contents = line.split(",")
        
        corrs = contents[1:numb_datasets+1]
        corrs_flag = 0
        
        if counter == 1: # move to next iteration
            of.write(line)
            continue
        
        # write some condition to check for NA
        pos = sum(float(num) > 0 for num in corrs)
        neg = sum(float(num) < 0 for num in corrs)
        
        #print pos, neg
        if len(corrs) == pos and not len(corrs) == neg:
            print "pos"
            corrs_flag = 1
        if len(corrs) == neg and not len(corrs) == pos:
            print "neg"
            corrs_flag = 1

        if corrs_flag == 1: 
            corrs_dict[contents[0]] = contents[1:]
        '''
        if corrs_flag == 0: # no point in analyzing pvals, move to next iteration
            continue

        pvals = contents[numb_datasets+1:]
        print pvals
        pvals_flag = 0
                
        # write some condition to check for NA
        sig = sum(float(num) < 1 for num in pvals)
        
        #print sig
        if len(corrs) == sig:
            print "sig"
            pvals_flag = 1

        if corrs_flag == 1 and pvals_flag == 1:
            corrs_dict[contents[0]] = contents[1:]
        
        if counter == 5:
            sys.exit(0)
        '''
        

print corrs_dict
"""