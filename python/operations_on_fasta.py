import os
import sys
from utils import *
import subprocess
import argparse
from time import localtime, strftime

# usage: python ~/scripts/python/operations_on_fasta.py -i chimeras.fa -o chimeraSeqlength.txt -c
# python ~/scripts/python/operations_on_fasta.py -i derep_otus97.fa -o derep_otus97_seqs_filtered.fa -n 220 -x 280 -f

def main(args):

    parser = argparse.ArgumentParser(description='Calculate stats. on fasta file and some basic operations such as filtering seqs. of certain length.')
    parser.add_argument('-i', '--infile')
    parser.add_argument('-o', '--outfile', default="topk_taxon_out.txt")  # output filename
    parser.add_argument('-n', '--minlength', default=10, type=int) # keep seqs with length more than equal to minlength
    parser.add_argument('-x', '--maxlength', default=1000, type=int) # keep seqs with length less than equal to maxlength
    parser.add_argument('-c', '--calcstats', action='store_true', default=False) # calc. stats on the fasta file
    parser.add_argument('-f', '--filterseqs', action='store_true', default=False) # filter. seqs. that have seq length < minlength or > maxlength
    parser.add_argument('-d', '--delimiter', default='\t')  # delimiter for file
    
    args = parser.parse_args()

    if len(sys.argv)==1 :
        parser.print_help()
        sys.exit('\natleast one argument required\n')

    infile  = args.infile
    outfile = args.outfile
    minl = args.minlength
    maxl = args.maxlength
    delim = args.delimiter
    
    if args.calcstats:
        calc_stats_fasta(infile, outfile)
    if args.filterseqs:
        filter_seqs(infile, outfile, minl, maxl)

def filter_seqs(infile, outfile, minl, maxl):
    header_list, seq_list, length_list = create_header_seq_and_length_lists(infile)
    f = open(outfile, 'w')
    for idx, seq in enumerate(seq_list):
        if len(seq) >= minl and len(seq) <= maxl:
            f.write('%s\n' %header_list[idx])
            # split the seqs in 80 chars length and then write to file
            chunk_size=80
            seq_chunk_size_list = [ seq[i:i+chunk_size] for i in range(0, len(seq), chunk_size) ]
            f.write('%s\n' %('\n'.join(seq_chunk_size_list)))
    f.close()


def calc_stats_fasta(infile, outfile): # calc stats on the fasta file
    header_list, seq_list, length_list = create_header_seq_and_length_lists(infile)
    writeLIST_to_file(length_list, outfile)
    cmd='Rscript /home/williamslab/scripts/R/plothist.R ' + outfile
    subprocess.call(cmd.split())


def create_header_seq_and_length_lists(infile):
    seq_list = list()
    length_list = list()
    header_list = list()
    string = ''
    for l in open(infile):
        if '#' in l:
            continue
        else:
            l = l.strip()
        if '>' in l:
            header_list.append(l)
            if string != '':
                seq_list.append(string)
                length_list.append(len(string))
                string = ''
        else:
            string += l
    seq_list.append(string)
    length_list.append(len(string))
    #del seq_list[0]  # to remove the blank string added at the first time
    #del length_list[0]
    #print len(header_list) , len(seq_list), len(length_list)
    return header_list, seq_list, length_list


if __name__=='__main__':
    datetime = strftime("%a, %d %b %Y %I:%M:%S %p", localtime())
    cmd = 'echo ' + datetime 
    os.system(cmd)
    main(sys.argv)

