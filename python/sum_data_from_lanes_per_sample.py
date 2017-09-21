import os
import sys
from utils import *
import operator
from time import localtime, strftime
import argparse
import re
import sh
import pandas as pd


# cd /nfs3/PHARM/Morgun_Lab/richrr/Type2_Diabetes/RNA-Seq/analysis/summarize_per_sample
# usage:  python ~/Morgun_Lab/richrr/scripts/python/sum_data_from_lanes_per_sample.py -i data_per_lane.txt -b rnaseqil


# sort it first on the prefix and then suffix of the "strng"
def sort_headers(names, strng):
    return sorted(names, key=lambda x: (int(x.split(strng)[0]), int(x.split(strng)[1]))) 


def main(args):
    # extract the statistics from the log
    parser = argparse.ArgumentParser(description='sum the number of reads from different lanes per sample. The current delimiter is _L00')
    parser.add_argument('-i', '--infile')
    parser.add_argument('-o', '--outfile', default="summed-values.txt")  # output filename
    parser.add_argument('-d', '--delimiter', default='\t')  # delimiter for file
    parser.add_argument('-l', '--lanedelimit', default='_L00')  # delimiter for identifying lanes
    parser.add_argument('-r', '--rowid', default=0)  # row id
    parser.add_argument('-b', '--buffr')  # the buffer string is used to split the sample names and used
                                           #  during arranging samples order in output file
                                           # for samples named 11rnaseq1, 11rnaseq2, 12rnaseq1, ...the -b is rnaseq
                                           # the samples are first sorted using prefix of rnaseq and then suffix of rnaseq

    args = parser.parse_args()

    if len(sys.argv)==1 :
        parser.print_help()
        sys.exit('\natleast one argument required\n')

    infile  = args.infile
    buffr = args.buffr
    outfile = args.outfile
    delim = args.delimiter
    lanedelimit = args.lanedelimit
    
    # get the header line from the big file
    header_line = ''
    with open(infile, 'r') as f:
        header_line = f.readline().strip()
    #print header_line

    # create a dict on the sample name and the different lanes to be summed
    sample_lanes_dict =  dict()
    for entry in header_line.split(delim)[1:]:  # here we are ignoring the header of the first column (e.g. ID, SampleID)
        entry = entry.strip()

        # split the entry with the lane delimiter and fill the dictionary
        contents = entry.split(lanedelimit)
        key = contents[0]
        #print entry
        if key in sample_lanes_dict:
            sample_lanes_dict[key].append(entry)
        else:
            sample_lanes_dict[key] = [entry]
    #print sample_lanes_dict
   
    # sort the samples using the prefix and suffix of the buffr string
    headers = sort_headers(sample_lanes_dict.keys(), buffr)

    # read the full file as a data frame
    df = pd.read_table(infile, index_col = args.rowid)
    print df.head()

    # create an empty df with the required row and column names
    df_ = pd.DataFrame(index=list(df.index), columns=headers)
    df_ = df_.fillna('NaN')
    #print df_.head()

    # sum the specific columns of the data frame
    # https://stackoverflow.com/questions/25748683/pandas-sum-dataframe-rows-for-given-columns
    for key in headers:
        value = sample_lanes_dict[key]
        print key, value
        df_[key] = df[value].sum(axis=1)
    print df_.head()

    df_.to_csv(outfile, sep=delim, index_label=header_line.split(delim)[0])

    # print the column names
    print "Check the order of samples" 
    print '\n'.join(list(df_))
    print "Done"
    #sys.exit()

    

if __name__=='__main__':
    datetime = strftime("%a, %d %b %Y %I:%M:%S %p", localtime())
    cmd = 'echo ' + datetime 
    os.system(cmd)
    main(sys.argv)

