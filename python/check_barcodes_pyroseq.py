import os
import sys
from utils import *
import operator
from time import localtime, strftime
import argparse
import re
import os.path
from subprocess import Popen, PIPE
from collections import defaultdict
import numpy as np
import matplotlib.pyplot as plt
import math
#from Bio.Seq import Seq
import Bio.Seq as BioSeq


#usage:
#python ~/scripts/python/check_barcodes_pyroseq.py -i ../../analysis/michigan-illumina.txt -o barcode_stats_q30.csv -q 30 -z ./Trimmed_fastq_files/q30/

#python ~/scripts/python/check_barcodes_pyroseq.py -i ../../analysis/michigan-illumina.txt -o barcode_stats_q30.csv -q 30 -z ./Trimmed_fastq_files/q30/ -c


'''
Currently modified to work with Jesus et al. GCB bioenergy 2015 data
'''

def draw_histogram_adpaters_trimed_seqs(infile):
    lines = read_file(infile)
    counter = flag = allow_error_flag = leng_freq_flag = 0
    adapter_trimmed_dict = defaultdict(list)  # keeps track of how many times the adapter was cut without concern of length of reads ----(1)
    allowed_error_dict = dict()
    leng_freq_dict = defaultdict(int) # keeps track of how many reads were cut along with their size ------(2)
    leng_freq_per_adap_dict = defaultdict(int)
    adaptername = ''
    for l in lines:
        l = l.strip()
        #print l
        if not l:
            counter = allow_error_flag = leng_freq_flag = 0
            continue
        counter += 1
        if re.search('Adapter.*length.*was trimmed.*times', l) != None:
            ad, adaptername_, seq, l, len, w, t, freq, times = l.split(' ')   #----(1)
            adaptername = adaptername_.replace("'", "")
            #print l, adaptername, seq, len, freq
            adapter_trimmed_dict[adaptername].append(freq) #because of defaultdict(list), When each key is encountered for the first time, it is not already in the mapping; so an entry is automatically created using the default_factory function which returns an empty list.
            counter = 0
        if 'No. of allowed errors:' in l:
            counter = 0
            allow_error_flag = 1
        if allow_error_flag == 1 and counter == 1 and re.search('^0-19 bp: ', l) != None:
            #print l
            allowed_error_dict[adaptername] = l
            counter = 0
            allow_error_flag = 0
            #sys.exit(0)
        if 'length	count	expect	max.err	error counts' in l:
            counter = 0
            leng_freq_flag = 1
        if leng_freq_flag == 1  and counter >= 1:
            bases_removed, numb_reads, expec, max_err, err_counts = l.split('\t')  #----(2)
            #print bases_removed, numb_reads  # err_counts are space delimited, others are tab-delim
            leng_freq_dict[int(bases_removed)] += int(numb_reads)
            leng_freq_per_adap_dict[adaptername, int(bases_removed)] += int(numb_reads)  

    #print allowed_error_dict , '\n', leng_freq_dict , '\n' , leng_freq_per_adap_dict
    all_adaps_sum = numb_non_zero_values = numb_zero_values = avg = 0
    for k,v in adapter_trimmed_dict.items():
        l = [int(i) for i in v]   # convert entire list to int
        suml = sum(l)             # sum of list
        all_adaps_sum += suml 
        numb_non_zero_values = sum(x > 0 for x in l)   # number of non zero elements
        numb_zero_values = sum(1 for x in l if x == 0)      # number of zero elements
        #print k, l, suml, numb_non_zero_values, numb_zero_values
        avg = 0
        if numb_non_zero_values != 0:
            avg = suml/numb_non_zero_values
        text ='Adapter: %s trimmed an average of %d reads from %d samples, sample counted if atleast one of \n its read was trimmed for the adapter, from total %d samples,  allowed error %s' %(k, avg, numb_non_zero_values, numb_non_zero_values+numb_zero_values, allowed_error_dict[k])
        mydict = {x[1]:y for x,y in leng_freq_per_adap_dict.iteritems() if x[0]==k}
        #print text, '\n', mydict
        if mydict:
            for t in ['scatter', 'line']:
                draw_hist(t , k, mydict, k, text)
    avg = all_adaps_sum/numb_non_zero_values
    text ='All Adapter: trimmed an average of %d reads from %d samples, sample counted if atleast one of \n its read was trimmed for the adapter, from total %d samples,  allowed error ' %(avg, numb_non_zero_values, numb_non_zero_values+numb_zero_values)
    for a,e in allowed_error_dict.items():
        text += '%s %s \n' %(a,e)
    for t in ['scatter', 'line']:
        draw_hist(t , 'all_adapters', leng_freq_dict, 'all_adapters', text, -0.12)
    sys.exit(0)


def draw_hist(fig_type, fig_name, counted_data, fig_spec, text, ypos=-0.03):
    k, v = counted_data.keys(), counted_data.values() # the orders are unchanged
    if fig_type == 'line':
        plt.plot(k,v)
    if fig_type == 'scatter':
        plt.scatter(k,v)
    plt.yscale('log')
    # to draw histogram
    #plt.hist(counted_data.keys(), weights=counted_data.values(), bins=range(max(counted_data) + 10), log=True
    plt.title('# of trimmed bases vs. # of reads ' + fig_spec)
    plt.xlabel('Number of bases removed by trimming')
    plt.ylabel('Number of reads (log scale)')
    plt.figtext(0,ypos,text, fontsize='x-small')
    plt.savefig(fig_name+fig_type+'.png', bbox_inches='tight') #replace with .pdf
    plt.clf()

def grep_string(string, f):
    cmd = 'grep -c "^%s" %s' %(string,f) # pattern
    output_start = Popen(cmd, shell=True, stdout=PIPE).communicate()[0]  # grep -c 'pattern' filename
    #print cmd , '\n', output_start
    cmd = 'grep -c "%s$" %s' %(string,f) # pattern
    output_end = Popen(cmd, shell=True, stdout=PIPE).communicate()[0]  # grep -c 'pattern' filename
    #print cmd , '\n', output_end
    return output_start, output_end

def grep_string_contains(string, f):
    cmd = 'grep -c "%s" %s' %(string,f) # pattern
    output = Popen(cmd, shell=True, stdout=PIPE).communicate()[0]  # grep -c 'pattern' filename
    #print cmd , '\n', output
    return output

'''
def grep_two_patterns(start, end, f):
    cmd = 'egrep -c "^%s.+%s" %s' %(start, end,f) # pattern
    output = Popen(cmd, shell=True, stdout=PIPE).communicate()[0]  # grep -c 'pattern' filename
    #print cmd , '\n', output

    cmd = 'egrep -c "^%s.+%s$" %s' %(start, end,f) # pattern
    output = Popen(cmd, shell=True, stdout=PIPE).communicate()[0]  # grep -c 'pattern' filename
    print cmd , '\n', output

    return output

def grep_two_patterns_negate_first(start, end, f):
    cmd = 'grep -v "^%s" %s | grep -c "%s"' %(start, f, end) # pattern
    output = Popen(cmd, shell=True, stdout=PIPE).communicate()[0]  # grep -c 'pattern' filename
    #print cmd , '\n', output
    return output

def grep_two_patterns_negate_second(start, end, f):
    cmd = 'grep "^%s" %s | grep -v -c "%s"' %(start, f, end) # pattern
    output = Popen(cmd, shell=True, stdout=PIPE).communicate()[0]  # grep -c 'pattern' filename
    #print cmd , '\n', output
    return output

def grep_two_patterns_negate_both(start, end, f):
    cmd = 'grep -A 1 "^@M" %s | grep -v "^%s" | grep -v "^@" | grep -v "\-\-" | grep -v -c "%s"' %(f, start, end) # pattern
    output = Popen(cmd, shell=True, stdout=PIPE).communicate()[0]  # grep -c 'pattern' filename
    #print cmd , '\n', output
    return output
'''

def main(args):

    parser = argparse.ArgumentParser(description='Check if the barcodes  are present in the sample')
    parser.add_argument('-i', '--infile') # (sample info.) file containing the sample name and the barcode string
    parser.add_argument('-o', '--outfile', default="out.txt")  # output filename
    parser.add_argument('-l', '--logfile', default="ITS-log-file.txt")  # output filename
    parser.add_argument('-a', '--amplicon16s', action='store_true', default=False) # 16S, default ITS
    parser.add_argument('-f', '--fasta', action='store_true', default=False) # files in fasta format, default fastq
    parser.add_argument('-q', '--quality', default=30, type=int) # Phred quality threshold of seq.
    parser.add_argument('-z', '--trimfastq', default='./Trimmed_fastq_files/')  # prefix for trimmed files
    parser.add_argument('-c', '--cutadapt', action='store_true', default=False) # run cutadapt
    parser.add_argument('-d', '--delimiter', default='\t')  # delimiter for file
    parser.add_argument('-m', '--histogram', action='store_true', default=False) # draw histogram for each adapter using infile

    args = parser.parse_args()

    if len(sys.argv)==1 :
        parser.print_help()
        sys.exit('\natleast one argument required\n')

    infile  = args.infile
    outpfile = args.outfile
    filetype = '.fasta' if args.fasta else '.fastq'
    logfile = args.logfile
    run_cutadapt = args.cutadapt
    delim = args.delimiter
    if args.histogram:
        draw_histogram_adpaters_trimed_seqs(infile)
        sys.exit(0)

    global FWDPRIMER, FADAPTER, REVPRIMER, RADAPTER, QUALITY, R_1, R_2, INTERMFASTQ, INTERMINFO, TRIMFASTQF, TRIMFASTQD, \
                REV_COMPLEM_RAD_RP, REV_COMPLEM_FAD_FP, FAD_FP, RAD_RP, REV_COMPLEM_FWDPRIMER, REV_COMPLEM_REVPRIMER
    if args.amplicon16s: # primers and adapters for the 16S data
        FWDPRIMER='AAACTYAAAKGAATTGRCGG'
        FADAPTER='CCTATCCCCTGTGTGCCTTGGCAGTCTCAG'       
        REVPRIMER='ACGGGCGGTGTGTRC'
        RADAPTER='CCATCTCATCCCTGCGTGTCTCCGACTCAG'
    else: # primers and adapters for the ITS data
        pass
	#FWDPRIMER='CTTGGTCATTTAGAGGAAGTAA'
        #FADAPTER='TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG'       
        #REVPRIMER='GCTGCGTTCTTCATCGATGC'
        #RADAPTER='GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG'

    REV_COMPLEM_FWDPRIMER = BioSeq.reverse_complement(FWDPRIMER) # search at end in read 2
    REV_COMPLEM_REVPRIMER = BioSeq.reverse_complement(REVPRIMER) # search at end in read 1

    FAD_FP = FADAPTER + FWDPRIMER # check at start in read 1
    RAD_RP = RADAPTER + REVPRIMER # check at start in read 2
    REV_COMPLEM_RAD_RP = BioSeq.reverse_complement(RAD_RP) # search at end in read 1
    REV_COMPLEM_FAD_FP = BioSeq.reverse_complement(FAD_FP) # search at end in read 2

    #print FADAPTER, FWDPRIMER, FAD_FP, REV_COMPLEM_FAD_FP, REVPRIMER, RADAPTER, RAD_RP, REV_COMPLEM_RAD_RP, REV_COMPLEM_FWDPRIMER, REV_COMPLEM_REVPRIMER

    QUALITY=args.quality
    TRIMFASTQD = args.trimfastq  # dir
    TRIMFASTQF = TRIMFASTQD + 'trim_' # dir + prefix for files
    
    for dir in [TRIMFASTQD]:
        if not os.path.exists(dir):
            os.makedirs(dir)
    
    lines = read_file(infile)
    SAMPLE_BARCODE_DICT = dict()  # this is not a global variable
    sampleid_col = fbarcode_col = rbarcode_col = 0
    for l in lines:
        cont = l.strip().split(delim)
        if '#' in l:
            # identify column number of SampleId and barcodeseq. check if present to avoid ValueError if the target is not found
            if "#SampleID" in cont:
                sampleid_col = cont.index("#SampleID")
            continue
        if sampleid_col == -1:
            sys.exit('Could not detect the column labels SampleID, FBarcodeSequence, or RBarcodeSequence')
        # adding the suffix helps avoid processing the file twice since empty line/string
        # in samplesheet can still satisfy "if new_f.startswith(this_dir+trimk)"
        sample = cont[sampleid_col] + filetype
        if sample in SAMPLE_BARCODE_DICT:
            sys.exit('\nduplicate samples detected, check sample sheet\n')
        else:
            SAMPLE_BARCODE_DICT[sample] = 1
        # The pyroseq 16S metagenomic kit has the following layout:
        #fwd_adapter - fwd_primer - DNA - rev_primer - barcode - rev_adapter
    barcode_stats(outpfile, SAMPLE_BARCODE_DICT, delim, filetype, run_cutadapt, logfile)

'''
R = G A (purine)
Y = T C (pyrimidine)
K = G T (keto)
M = A C (amino)
S = G C (strong bonds)
W = A T (weak bonds)
B = G T C (all but A)
D = G A T (all but C)
H = A C T (all but G)
V = G C A (all but T)
N = A G C T (any)
badset = set(R,Y,K,M,S,W,B,D,H,V,N)
'''

#http://stackoverflow.com/questions/20726010/how-to-check-if-a-string-contains-only-characters-from-a-given-set-in-python
def replace_bad_nucleotides(pattern):
    good_set = set('ATGC')
    bad_set = set('RYKMSWBDHVN')
    if set(pattern) <= good_set:
        pass
    else:
        for i in pattern:
            if i in bad_set:
                pattern = pattern.replace(i, '.')
    #print pattern
    return pattern


def barcode_stats(outpfile, SAMPLE_BARCODE_DICT, delim, filetype, run_cutadapt, logfile, this_dir='./', trimmedbarcodestatus = False):
    # you do not need to mention global to read values of global variables, only if you want to re/assign value
    outfile = open(outpfile, 'w')
    headers_list = ['file' ,\
    'starts_FADAPTER', 'ends_FADAPTER','contains_FADAPTER','starts_FWDPRIMER', 'ends_FWDPRIMER', 'contains_FWDPRIMER',\
    'starts_FAD_FP', 'ends_FAD_FP', 'contains_FAD_FP',\
    'starts_REV_COMPLEM_FAD_FP', 'ends_REV_COMPLEM_FAD_FP', 'contains_REV_COMPLEM_FAD_FP',\
    'starts_REVPRIMER', 'ends_REVPRIMER', 'contains_REVPRIMER',\
    'starts_RADAPTER', 'ends_RADAPTER', 'contains_RADAPTER','starts_RAD_RP', 'ends_RAD_RP', 'contains_RAD_RP',\
    'starts_REV_COMPLEM_RAD_RP', 'ends_REV_COMPLEM_RAD_RP', 'contains_REV_COMPLEM_RAD_RP',\
    'starts_REV_COMPLEM_FWDPRIMER', 'ends_REV_COMPLEM_FWDPRIMER', 'contains_REV_COMPLEM_FWDPRIMER',\
    'starts_REV_COMPLEM_REVPRIMER', 'ends_REV_COMPLEM_REVPRIMER', 'contains_REV_COMPLEM_REVPRIMER']
    outfile.write ("%s\n" % delim.join(headers_list))
    for f in sorted(os.listdir(this_dir)):
        new_f = this_dir + f
        if os.path.isfile(new_f) and new_f.endswith(filetype):
            for k in SAMPLE_BARCODE_DICT:  # for each element in dict
                trimk = k
                if (trimmedbarcodestatus):
                    trimk = 'trim_' + k
                # uncomment these next two lines if you directly want to run the barcode status on trimmed reads
                #else:
                #    continue                
                if new_f.startswith(this_dir+trimk): #if this_dir+trimk in new_f:   # if the key is prefix in filename.
                    #print 'ok', '--', trimk, '--' , new_f
                    counts_list = [f]
                    for pattern in [FADAPTER, FWDPRIMER, FAD_FP, REV_COMPLEM_FAD_FP,\
                                      REVPRIMER, RADAPTER, RAD_RP, REV_COMPLEM_RAD_RP, REV_COMPLEM_FWDPRIMER, REV_COMPLEM_REVPRIMER]:
                        pattern = replace_bad_nucleotides(pattern)
                        counts_list.extend(grep_string(pattern, new_f)) # separately checks starts with pattern, ends with pattern
                        counts_list.append(grep_string_contains(pattern, new_f))

                    str_counts_list = [str(i.strip()) for i in counts_list]
                    outfile.write ("%s\n" % delim.join(str_counts_list))
    outfile.close()
    if (run_cutadapt):
        print "Done with initial barcode stats\nNow running cutadapt\n"
        run_cutadapt_method(outpfile, logfile, this_dir, filetype, SAMPLE_BARCODE_DICT, delim)
    else:
        sys.exit('\nRunning cutadapt was not requested.\n')


def run_cutadapt_method(outpfile, logfile, this_dir, filetype, SAMPLE_BARCODE_DICT, delim):
    completed_files_list = list() 
    outfile = open(logfile, 'w')
    for f in sorted(os.listdir(this_dir)):
        if os.path.isfile(f) and f.endswith(filetype):
            for k in SAMPLE_BARCODE_DICT:  # for each element in dict
                if f == k and f not in completed_files_list:                         ## fix this condition
                    # if the key is prefix in filename. did this instead of k in f 
                    #to avoid 2_S being true 72_S, 62_S, etc. and causing duplicates
                    #put the read 1 and read 2 file into completed list
                    #go ahead and do the trimming
                    completed_files_list.append(f)
                    # trim, qf
                    output_dump = quality_trim_min_length (f)
                    outfile.write ("%s\n---- Sample ----\n" % ('\n'.join(output_dump)))
    outfile.close()
    outpfile = "trimmed_" + outpfile
    # set run_cutadapt arg to false and run barcode stats
    print "\nDone with cutadapt\nNow running final barcode stats\n"
    run_cutadapt = False
    barcode_stats(outpfile, SAMPLE_BARCODE_DICT, delim, filetype, run_cutadapt, logfile, this_dir=TRIMFASTQD, trimmedbarcodestatus=True)


def quality_trim_min_length (f):
    outfile = TRIMFASTQF + f
    infofile = TRIMFASTQF + f + '_info.txt'
    min_length = 150
    RIDIC_ADAPTER = 'ADAPTER' # to get around the requirement of providing a true adapter; needs to be given with '-N'
    #suffix_string = 'orig_bc=XXXX_new_bc=XXXX_bc_diffs=0' # gives "too many parameters error" if the string has spaces
    #cmd = 'cutadapt -q %s --minimum-length %s --overlap=10 --info-file=%s -N -a %s --suffix %s -o %s %s' %(QUALITY, min_length, infofile, RIDIC_ADAPTER, suffix_string, outfile, f)
    cmd = 'cutadapt -q %s --minimum-length %s --overlap=10 --info-file=%s -N -a %s -o %s %s' %(QUALITY, min_length, infofile, RIDIC_ADAPTER, outfile, f)
    output = Popen(cmd, shell=True, stdout=PIPE).communicate()[0]
    #print cmd, '\n', output
    qiime_compatible(f) # convert to qiimeformat and fasta
    return output


def qiime_compatible(f):
    sname = f.replace('.fastq','')
    outsname = TRIMFASTQF + f
    outputfile = TRIMFASTQF + sname
    barcode='XXXX'
    cmd='awk -f ~/scripts/awk/qiime_redo_headers_from_fastq_to_fasta.awk ">%s" "orig_bc=%s new_bc=%s bc_diffs=0" %s >%s.fasta' \
                                                      %(sname, barcode, barcode, outsname, outputfile)
    #print cmd
    #sys.exit(0)
    output_start = Popen(cmd, shell=True, stdout=PIPE).communicate()[0]


if __name__=='__main__':
    datetime = strftime("%a, %d %b %Y %I:%M:%S %p", localtime())
    cmd = 'echo ' + datetime 
    os.system(cmd)
    main(sys.argv)

