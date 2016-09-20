import os
import sys
from utils import *
import operator
from time import localtime, strftime
import argparse
import re
from Bio import SeqIO
from abifpy import Trace
import pylab
# using abifpy because it allows to change the cutoff value for quality trimming

# Usage: python ~/scripts/python/sanger_seq_analysis.py -i ../ -o se_v0_july -s se_v0_july
#        python ~/scripts/python/sanger_seq_analysis.py -i ../ -o se_v0_july -s se_v0_july -e fasta

def main(args):
    # read sanger seq files
    # output to format
    # blast
    parser = argparse.ArgumentParser(description='Analyze Sanger seq abi files')
    parser.add_argument('-i', '--indir', default="./") # directory containing ab1 files to be processed
    parser.add_argument('-o', '--outfile', default="default_output_file")  # output filename
    #parser.add_argument('-g', '--groupsfile')  # file containing the group info for each sample
    #parser.add_argument('-c', '--categories' , nargs='*')  \
              # creates a list of groups (categories) from the group file to calc. stats for, currently only one category
    parser.add_argument('-d', '--delimiter', default='\t')  # delimiter for file
    parser.add_argument('-e', '--extension', default='fastq')  # output file format (extension) options allowed: fasta
    parser.add_argument('-l', '--minseqlength', default=250, type=int) # minimum sequence length before trimming
    parser.add_argument('-q', '--phredqthreshold', default=25, type=int) # phred quality threshold, 1 error in approx. 316 bases
    # The fastx toolkit quality filter and trim method is a different strategy
    # I am using the quality trimming strategy as used by afifpy, clc and biopython.
    # Trims the sequence using Richard Mott's modified trimming algorithm.
    # Note: cutadapt has a similar strategy (which I used for the quality filtering the read 2 from fungal invasive project)
    # http://www.clcsupport.com/clcgenomicsworkbench/650/Quality_trimming.html
    # http://www.sourcecodebrowser.com/python-biopython/1.59/namespace_bio_1_1_seq_i_o_1_1_abi_i_o.html#a7fb9b10aa5e88c3f1da16b212f364a3a
    parser.add_argument('-r', '--revprimer', default='1492R')  # name of reverse primer
    parser.add_argument('-s', '--fillerstring', default='SWG__') # string to make the header name >10 chars. long to make libshuff compatible
    
    args = parser.parse_args()

    print "\nPlease double check if correct name of reverse primer is in the ab1 file names\n"

    if len(sys.argv)==1 :
        parser.print_help()
        sys.exit('\natleast one argument required\n')

    indir  = args.indir
    outfile = args.outfile
    delim = args.delimiter
    extens = args.extension
    segment = args.minseqlength

    phred_q_threshold = args.phredqthreshold
    cutoff = 10 ** (phred_q_threshold / -10.0)    # base-calling error probabilities as the trim threshold
    # Bio.SeqIO.AbiIO._abi_trim only uses default phred score 13, which is roughly 10^-1.3 = 0.05 base-calling error probabilities as the trim threshold
    # so copied the method and added the cutoff parameter in the code
    # so I can use phred score 25, i.e. 10^-2.5 ~ 0.0032 base-calling error probabilities as the trim threshold
    # since phred score 30, 10^-3 might be too strict

    rev_primer = args.revprimer
    filler_string = '>' + args.fillerstring + '_'


    #read_abif_file(indir,outfile, delim, extens, cutoff, phred_q_threshold)
    infile_for_groups_file = read_abif_file_biopy(indir,outfile, delim, extens, cutoff, phred_q_threshold, segment, rev_primer, filler_string)
    create_groups_file(infile_for_groups_file, filler_string, extens)


def read_abif_file_biopy(indir,outfile, delim, extens, cutoff, phred_q_threshold, segment, rev_primer, filler_string):
    outputfile_orig = outfile + "." + extens
    outputfile_trim = outfile + "_qtrim_phred_" + str(phred_q_threshold) + "." + extens
    outputfile_fwd_dir = outfile + "_fwd_dir_qtrim_phred_" + str(phred_q_threshold) + "." + extens
    
    if os.path.isfile(outputfile_orig) or os.path.isfile(outputfile_trim) or os.path.isfile(outputfile_fwd_dir):
        sys.exit("\nOutput file(s) already exist(s)! Did not run further code.\n")
    fh_orig = open(outputfile_orig, "a")
    fh_trim = open(outputfile_trim, "a")
    fh_fwd_dir = open(outputfile_fwd_dir, "a")

    ls = os.listdir(indir)
    for file in ls:
        if file[-4:] == '.ab1':
            for seq_record in SeqIO.parse(indir+file, "abi"):
                fh_orig.write(seq_record.format(extens))

                # trim the seqs based on phred quality threshold
                qual_trimmed_seq_record = richrr_hack_abi_trim(seq_record, cutoff, segment)
                fh_trim.write(qual_trimmed_seq_record.format(extens))

                rec_to_print = ''
                # another file, reverse complement the reverse primer seqs after quality filtering
                if rev_primer in qual_trimmed_seq_record.name:
                    rev_complem_qual_trimmed_seq_record = qual_trimmed_seq_record.reverse_complement(id="RC_"+qual_trimmed_seq_record.id, description = "reverse_complement")
                    #print qual_trimmed_seq_record.letter_annotations, '\n\n', rev_complem_qual_trimmed_seq_record.letter_annotations
                    rec_to_print = rev_complem_qual_trimmed_seq_record.format(extens)
                else:
                    rec_to_print = qual_trimmed_seq_record.format(extens)
                
                if extens == 'fasta': # since I run libshuff using fasta, this is only for fasta
                    rec_to_print = rec_to_print.replace('>',filler_string)
                    rec_to_print = rec_to_print.replace(' ','_')
                fh_fwd_dir.write(rec_to_print)

                #print seq_record.seq, '\n\n', qual_trimmed_seq_record.seq

    fh_orig.close()
    fh_trim.close()
    return outputfile_fwd_dir

def create_groups_file(infile, filler_string, extens):
    if extens == 'fastq':
        print_qual_scores_plot(infile)
    if extens == 'fasta':
        print_seq_length_stats_plot(infile)
        groups_file = 'groups-file.txt'
        header_groups_list = list()
        lines = read_file(infile)  # file containing fasta seqs.
        for l in lines:
            l = l.strip()
            if '>' in l:
                g = ''
                if 'KXA' in l: #'_A' 'KA' 'KXA'
                    g = 'A'
                elif 'KXD' in l: #'_D' 'KD' 'KXD'
                    g = 'D'
                str = l[1:] + '\t' + g
                header_groups_list.append(str)
        writeLIST_to_file(header_groups_list, groups_file)
        print '\nCheck the header and groups in %s\n' %groups_file

    print '\nCan only create groups file and seq length plot for fasta format\n \
        Can only create quality score plot for fastq format\n'


def print_seq_length_stats_plot(fa_infile):
    sizes = [len(rec) for rec in SeqIO.parse(fa_infile, "fasta")]
    print len(sizes), min(sizes), max(sizes)
    pylab.hist(sizes, bins=20)
    pylab.title("%i sanger sequences\nLengths %i to %i" % (len(sizes),min(sizes),max(sizes)))
    pylab.xlabel("Sequence length (bp)")
    pylab.ylabel("Count")
    pylab.savefig(fa_infile+"_seq_len.pdf")


def print_qual_scores_plot(fq_infile):
    for i,record in enumerate(SeqIO.parse(fq_infile, "fastq")):
        if i >= 50 : break #trick!
        pylab.plot(record.letter_annotations["phred_quality"])
    pylab.ylim(0,45)
    pylab.ylabel("PHRED quality score")
    pylab.xlabel("Position")
    pylab.savefig(fq_infile+"_qual_score.pdf")



# the method below works, but if you need to write the original records, you need to add those
# many .write statements as below
# def read_abif_file(indir, outfile, delim, extens, cutoff, phred_q_threshold):
#    outputfile = outfile + "_qtrim_phred_" + str(phred_q_threshold) + "." + extens
#    if os.path.isfile(outputfile):
#        sys.exit("\nOutput file already exists! Did not run further code.\n")
#     fh = open(outputfile, "a")
#     ls = os.listdir(indir)
#     for file in ls:
#         if file[-4:] == '.ab1':
#             record = Trace(indir+file)
#             #record.export(outputfile, extens) # overwrites file
#             # trim outputs only seqs, unlike biopython which returns the trimmed seq record
#             # cannot read and trim with threshold at same time, Trace(indir+file, trimming=True, cutoff=0.001) is not allowed
#             qual_trimmed_record_seq = record.trim(record.seq, cutoff)
#             qual_trimmed_record_qual = record.trim(record.qual, cutoff)
#             if extens == 'fasta':
#                 #fh.write(record.seq)
#                 header = '>' + record.name + '\n'
#                 fh.write(header)
#                 fh.write(qual_trimmed_record_seq)
#                 fh.write('\n')
#             elif extens == 'fastq':
#                 header = '@' + record.name + '\n'
#                 fh.write(header)
#                 fh.write(qual_trimmed_record_seq)
#                 fh.write('\n+\n')
#                 fh.write(qual_trimmed_record_qual)
#                 fh.write('\n')    
#             sys.exit(0)
#     fh.close()


def richrr_hack_abi_trim(seq_record, cutoff, segment): 
    """ Added the cutoff and segment parameters so I can give values other than the default 0.05 and 20 resp
    Hacked http://biopython.org/DIST/docs/api/Bio.SeqIO.AbiIO-pysrc.html#_abi_trim
    Trims the sequence using Richard Mott's modified trimming algorithm.
 
    Arguments: 
        - seq_record - SeqRecord object to be trimmed.
        - cutoff
        - segment
 
    Trimmed bases are determined from their segment score, which is a 
    cumulative sum of each base's score. Base scores are calculated from 
    their quality values. 
 
    More about the trimming algorithm: 
    http://www.phrap.org/phredphrap/phred.html 
    http://www.clcbio.com/manual/genomics/Quality_abif_trimming.html 
    """ 
 
    start = False   # flag for starting position of trimmed sequence 
    trim_start = 0  # init start index 
 
    if len(seq_record) <= segment: 
        return seq_record 
    else: 
        # calculate base score 
        score_list = [cutoff - (10 ** (qual / -10.0)) for qual in seq_record.letter_annotations['phred_quality']] 
 
        # calculate cummulative score 
        # if cummulative value < 0, set it to 0 
        # first value is set to 0, because of the assumption that 
        # the first base will always be trimmed out 
        cummul_score = [0] 
        for i in range(1, len(score_list)): 
            score = cummul_score[-1] + score_list[i] 
            if score < 0: 
                cummul_score.append(0) 
            else: 
                cummul_score.append(score) 
                if not start: 
                    # trim_start = value when cummulative score is first > 0 
                    trim_start = i 
                    start = True 
 
        # trim_finish = index of highest cummulative score, 
        # marking the end of sequence segment with highest cummulative score 
        trim_finish = cummul_score.index(max(cummul_score)) 
 
        return seq_record[trim_start:trim_finish] 
                    



if __name__=='__main__':
    datetime = strftime("%a, %d %b %Y %I:%M:%S %p", localtime())
    cmd = 'echo ' + datetime 
    os.system(cmd)
    main(sys.argv)

