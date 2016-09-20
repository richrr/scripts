import os
import sys
from utils import *
import operator
from time import localtime, strftime
import argparse
import re
from abifpy import Trace
from os import listdir
from os.path import isfile, join

# Usage: python ~/scripts/python/parseab1.py -i ./
# python ~/scripts/python/parseab1.py -f ./Well_A01_13R+1492R.ab1

def main(args):
    #https://pypi.python.org/pypi/abifpy/0.9
    parser = argparse.ArgumentParser(description='extracts sequence and various other data from Applied Biosystems, Inc. format (ABI) file.')
    parser.add_argument('-i', '--indir') # process all ab1 files in the directory
    parser.add_argument('-o', '--outfile', default="sanger_entire_seq_out.fa")  # output filename
    parser.add_argument('-f', '--infiles' , nargs='*')  # creates a list of files to be processed
    parser.add_argument('-d', '--delimiter', default='\t')  # delimiter for file
    
    args = parser.parse_args()

    if len(sys.argv)==1 :
        parser.print_help()
        sys.exit('\natleast one argument required\n')

    infiles = list()

    if args.indir != None:
        mypath = args.indir
        infiles = [ (join(mypath,f)) for f in listdir(mypath) if (isfile(join(mypath,f)) and 'ab1' in f) ]
    elif args.infiles != None:
        infiles =  args.infiles
    else:
        parser.print_help()
        sys.exit('\natleast one input required\n')
   
    outfile = args.outfile
    delim = args.delimiter

    out = open(outfile, 'w')
    for f in infiles:
        ab1file = Trace(f)
        out.write('>Sample-%s:File-%s\n%s\n' % (ab1file.name, ab1file.id, ab1file.seq))
        #ab1file.export(outfile)
    out.close()
        
         

if __name__=='__main__':
    datetime = strftime("%a, %d %b %Y %I:%M:%S %p", localtime())
    cmd = 'echo ' + datetime 
    os.system(cmd)
    main(sys.argv)

