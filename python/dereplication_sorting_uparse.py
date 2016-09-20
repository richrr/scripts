import os
import sys
from utils import *
import operator
from time import localtime, strftime
import argparse


# usage:
#python ~/scripts/python/dereplication_sorting_uparse.py -i ../test_sample_data_qiime_format.fasta
#python ~/scripts/python/dereplication_sorting_uparse.py -i ../test_sample_data_qiime_format.fasta -a
'''
williamslab@williamslab-vostro:~/Data/Ping$ ~/bin/usearch61 -derep_fulllength test_sample_data_qiime_format.fasta -output derep.fa -sizeout
usearch_i86linux32 v6.1.544, 4.0Gb RAM (24.7Gb total), 8 cores
(C) Copyright 2010-12 Robert C. Edgar, all rights reserved.
http://drive5.com/usearch

Licensed to: richrr@vt.edu

00:00 2.0Mb Reading test_sample_data_qiime_format.fasta, 16Mb
00:00  18Mb 50000 (50.0k) seqs, min 165, avg 254, max 300nt
00:00  80Mb 39523 (39.5k) uniques, avg cluster 1.3, median 1, max 274
00:00  80Mb  100.0% Writing derep.fa

'''

def main(args):

    parser = argparse.ArgumentParser(description='Given a set of sequences, count frequency of each unique sequence and sort as per the frequency. \
                                                   Default: writes seq (in desc. order) seen atleast twice to outputfile and singletons to outputfile_singletons file')
    parser.add_argument('-i', '--infile')
    parser.add_argument('-o', '--outfile', default="derep.fa")  # output filename. If specified the singletons will be written to filename_singletons.fa
    parser.add_argument('-n', '--nosingletons', action='store_true', default=False) # do not write the singletons to file
    parser.add_argument('-m', '--minsize', default=2, type=int) # Sequences seen at least these times are written to output file
    parser.add_argument('-a', '--ascendsort', action='store_true', default=False) # print the seq to file as per the ascending order of freq.
    parser.add_argument('-d', '--delimiter', default='\t')  # delimiter for file

    args = parser.parse_args()

    if len(sys.argv)==1 :
        parser.print_help()
        sys.exit('\natleast one argument required\n')

    infile  = args.infile
    outpfile = args.outfile
    min_size = args.minsize
    no_singletons = args.nosingletons
    asc_sort = args.ascendsort
    singletonsfile = 'singletons.fa'
    if outpfile != 'derep.fa':
        singletonsfile = outpfile.replace('.fa', '_singletons.fa')
    delim = args.delimiter
    
    '''
>A10_R12__0 orig_bc=GTATGCGCTGTA new_bc=GTATGCGCTGTA bc_diffs=0
TACGGGGGGGGCAAGCGTTGTTCGGAATTACTGGGCGTAAAGGGCGCGTAGGCGGTTCGCTAAGTTGAATGTGAAAACTCTGGGCTTAACCCAGAGCCTGCATTCAATACTGGCGGGCTAGAGTTCTGGAGGGGATAGCGGAATTCCTGGTGTAGCGGTGAAATGCGTAGATATCAGGAGGAACACCGGTGGCGAAGGCGGCTATCTGGACAGACTCTGACGCTGAGGCGCGAAAGCTAGGGGAGCGANCAGG
    '''
    header = ''
    # dictionary, where seq is key, value is a list of (first) header and the freq of occurance of seq
    seq_header_freq_dict = dict()
    with open(infile) as f:
        for l in f:
            l = l.strip()
            if '#' in l:
                continue
            if l[0] == '>':
                header = l
            else:
                seq = l
                if seq in seq_header_freq_dict: # increment freq of occurance of seq
                     seq_header_freq_dict[seq][1] += 1
                else: # observed the seq for the first time
                    seq_header_freq_dict[seq] = [header, 1]
    f.close()
    #print len(seq_header_freq_dict)
    
    # desc. order the dict as per the freq., dict[seq] = [header, freq]
    # http://stackoverflow.com/questions/1217251/python-sorting-a-dictionary-of-lists
    # descending (default) sort as per the last column
    sorted_seq_header_freq_dict = ''
    if asc_sort:
        sorted_seq_header_freq_dict = sorted(seq_header_freq_dict.items(), key=lambda k:k[1][1])
    else:
        sorted_seq_header_freq_dict = sorted(seq_header_freq_dict.items(), key=lambda k:k[1][1], reverse=True)
    
    f = open(outpfile, 'w')
    s = open(singletonsfile, 'w')
    for seq , value in sorted_seq_header_freq_dict:
        if value[1] >= min_size:  # print seq with atleast minsize freq
            header = ';'.join(value[0].split(' '))
            frequency = 'size=' + str(value[1])
            header_line = ';'.join([header, frequency])
            f.write('%s\n%s\n' %(header_line, seq))    
        elif no_singletons == False: # write singletons (or less than min_size) to singletons.fa
            header = ';'.join(value[0].split(' '))
            frequency = 'size=' + str(value[1])
            header_line = ';'.join([header, frequency])
            s.write('%s\n%s\n' %(header_line, seq))        
    f.close()
    s.close()



    '''
    f = open(outpfile, 'w')
    for seq in seq_header_freq_dict:
        header = ';'.join(seq_header_freq_dict[seq][0].split(' '))
        frequency = 'size=' + str(seq_header_freq_dict[seq][1])
        header_line = ';'.join([header, frequency])
        f.write('%s\n%s\n' %(header_line, seq))    

    f.close()
    '''

if __name__=='__main__':
    datetime = strftime("%a, %d %b %Y %I:%M:%S %p", localtime())
    cmd = 'echo ' + datetime 
    os.system(cmd)
    main(sys.argv)

