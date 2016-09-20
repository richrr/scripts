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
#python ~/scripts/python/check_barcodes_nextera.py -i ../../analysis/michigan-illumina.txt -o barcode_stats_q30.csv -q 30 -z ./Trimmed_fastq_files/q30/

#python ~/scripts/python/check_barcodes_nextera.py -i ../../analysis/michigan-illumina.txt -o barcode_stats_q30.csv -q 30 -z ./Trimmed_fastq_files/q30/ -c

    

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
    parser.add_argument('-i', '--infile') # file containing the sample name and the barcode string
    parser.add_argument('-o', '--outfile', default="out.txt")  # output filename
    parser.add_argument('-l', '--logfile', default="ITS-log-file.txt")  # output filename
    parser.add_argument('-a', '--amplicon', default="") # 16S, nifH, ITS
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
    if args.amplicon == "16S": # primers and adapters for the 16S data
        FWDPRIMER='CCTACGGGNGGCWGCAG'
        FADAPTER='TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG'       
        REVPRIMER='GACTACHVGGGTATCTAATCC'
        RADAPTER='GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG'
    elif args.amplicon == "nifH": # primers and adapters for the nifH data
        FWDPRIMER='TGCGAYCCSAARGCBGACTC'
        FADAPTER='TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG'       
        REVPRIMER='ATSGCCATCATYTCRCCGGA'
        RADAPTER='GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG'
    elif args.amplicon == "ITS": # primers and adapters for the ITS data
        FWDPRIMER='CTTGGTCATTTAGAGGAAGTAA'
        FADAPTER='TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG'       
        REVPRIMER='GCTGCGTTCTTCATCGATGC'
        RADAPTER='GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG'
    else:
        sys.err("The code was modified so that you now have to explicitly mention the amplicon (16S, nifH, ITS). Previous code allowed 16S with a boolean (True) -a and without the argument it defaulted to ITS")

    REV_COMPLEM_FWDPRIMER = BioSeq.reverse_complement(FWDPRIMER) # search at end in read 2
    REV_COMPLEM_REVPRIMER = BioSeq.reverse_complement(REVPRIMER) # search at end in read 1

    FAD_FP = FADAPTER + FWDPRIMER # check at start in read 1
    RAD_RP = RADAPTER + REVPRIMER # check at start in read 2
    REV_COMPLEM_RAD_RP = BioSeq.reverse_complement(RAD_RP) # search at end in read 1
    REV_COMPLEM_FAD_FP = BioSeq.reverse_complement(FAD_FP) # search at end in read 2


    QUALITY=args.quality
    R_1 = '_L001_R1_001'
    R_2 = '_L001_R2_001'
    #INTERMFASTQ = args.intermedfastq
    #INTERMINFO = args.intermedinfo
    TRIMFASTQD = args.trimfastq  # dir
    TRIMFASTQF = TRIMFASTQD + 'trim_' # dir + prefix for files
    
    for dir in [TRIMFASTQD]:
        if not os.path.exists(dir):
            os.makedirs(dir)
    
    lines = read_file(infile)
    #GRL4463_S49_L001_R1_001 or 1_S72_L001_R1_001
    SAMPLE_BARCODE_DICT = dict()  # this is not a global variable
    sampleid_col = fbarcode_col = rbarcode_col = 0
    for l in lines:
        cont = l.strip().split(delim)
        if '#' in l:
            # identify column number of SampleId and barcodeseq. check if present to avoid ValueError if the target is not found
            if "#SampleID" in cont:
                sampleid_col = cont.index("#SampleID")
            if "FBarcodeSequence" in cont:
                fbarcode_col = cont.index("FBarcodeSequence")
            if "RBarcodeSequence" in cont:
                rbarcode_col = cont.index("RBarcodeSequence")
            continue
        if sampleid_col == 0 and fbarcode_col == 0 and rbarcode_col == 0:
            sys.exit('Could not detect the column labels SampleID, FBarcodeSequence, or RBarcodeSequence')
        sample = cont[sampleid_col] + '_S'
        barcode = cont[fbarcode_col] + '_' + cont[rbarcode_col]
        if sample in SAMPLE_BARCODE_DICT:
            sys.exit('\nduplicate samples/barcodes detected, check sample sheet\n')
        else:
            SAMPLE_BARCODE_DICT[sample] = barcode
        # The nextera 16S metagenomic kit has the following layout:
        #P5 - Index2 - overhang_adapter - fwd_primer - DNA - rev_primer - overhang_adapter - Index1 - P7
    barcode_stats(outpfile, SAMPLE_BARCODE_DICT, delim, filetype, run_cutadapt, logfile)


def barcode_stats(outpfile, SAMPLE_BARCODE_DICT, delim, filetype, run_cutadapt, logfile, this_dir='./', trimmedbarcodestatus = False):
    # you do not need to mention global to read values of global variables, only if you want to re/assign value
    outfile = open(outpfile, 'w')
    headers_list = ['file' , 'starts_fbarcode', 'ends_fbarcode', 'contains_fbarcode',\
    'starts_FADAPTER', 'ends_FADAPTER','contains_FADAPTER','starts_FWDPRIMER', 'ends_FWDPRIMER', 'contains_FWDPRIMER',\
    'starts_FAD_FP', 'ends_FAD_FP', 'contains_FAD_FP',\
    'starts_REV_COMPLEM_FAD_FP', 'ends_REV_COMPLEM_FAD_FP', 'contains_REV_COMPLEM_FAD_FP',\
    'starts_REVPRIMER', 'ends_REVPRIMER', 'contains_REVPRIMER',\
    'starts_RADAPTER', 'ends_RADAPTER', 'contains_RADAPTER','starts_RAD_RP', 'ends_RAD_RP', 'contains_RAD_RP',\
    'starts_REV_COMPLEM_RAD_RP', 'ends_REV_COMPLEM_RAD_RP', 'contains_REV_COMPLEM_RAD_RP',\
    'starts_rbarcode', 'ends_rbarcode' , 'contains_rbarcode',\
    'starts_REV_COMPLEM_fbarcode', 'ends_REV_COMPLEM_fbarcode' , 'contains_REV_COMPLEM_fbarcode',\
    'starts_REV_COMPLEM_rbarcode', 'ends_REV_COMPLEM_rbarcode' , 'contains_REV_COMPLEM_rbarcode',\
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
                    barcodes = SAMPLE_BARCODE_DICT[k].split("_")
                    #print 'ok', trimk, new_f
                    fbarcode = barcodes[0]
                    rbarcode = barcodes[1]
                    counts_list = [f]
                    rev_complem_fbarcode = BioSeq.reverse_complement(fbarcode) # search at end in read 2
                    rev_complem_rbarcode = BioSeq.reverse_complement(rbarcode) # search at end in read 1
                    for pattern in [fbarcode, FADAPTER, FWDPRIMER, FAD_FP, REV_COMPLEM_FAD_FP,\
                                      REVPRIMER, RADAPTER, RAD_RP, REV_COMPLEM_RAD_RP, rbarcode,\
                                      rev_complem_fbarcode, rev_complem_rbarcode, REV_COMPLEM_FWDPRIMER, REV_COMPLEM_REVPRIMER]:
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
                if f.startswith(k) and f not in completed_files_list:
                    # if the key is prefix in filename. did this instead of k in f 
                    #to avoid 2_S being true 72_S, 62_S, etc. and causinf duplicates
                    #put the read 1 and read 2 file into completed list
                    #go ahead and do the trimming
                    barcode = SAMPLE_BARCODE_DICT[k]
                    read1 = read2 = ''
                    if R_1 in f:
                        read1 = f
                        read2 = f.replace(R_1, R_2)
                    elif R_2 in f:
                        read2 = f
                        read1 = f.replace(R_2, R_1)
                    completed_files_list.extend([read1, read2])
                    # trim, qf
                    output_dump = remove_adap_bc_primer_quality_trim_min_length (read1, read2, barcode)
                    outfile.write ("%s\n---- Sample ----\n" % ('\n'.join(output_dump)))
    outfile.close()
    outpfile = "trimmed_" + outpfile
    # set run_cutadapt arg to false and run barcode stats
    print "\nDone with cutadapt\nNow running final barcode stats\n"
    run_cutadapt = False
    barcode_stats(outpfile, SAMPLE_BARCODE_DICT, delim, filetype, run_cutadapt, logfile, this_dir=TRIMFASTQD, trimmedbarcodestatus=True)


def remove_adap_bc_primer_quality_trim_min_length (read1, read2, barcode):
    #barcodes = barcode.split("_")
    #fbarcode = barcodes[0]
    #rbarcode = barcodes[1]
    outfile1 = TRIMFASTQF + read1
    outfile2 = TRIMFASTQF + read2
    infofile = TRIMFASTQF + read1 + "_" + read2 + '_info.txt'
    overlap = 10
    min_length = 100
    #The -q (or --trim-qualities) parameter can be used to trim low-quality ends from reads before adapter removal. 
    cmd = 'cutadapt -q %s --minimum-length %s --overlap %s --info-file=%s -g %s -G %s -a %s -A %s -n 3 -o %s -p %s %s %s' \
		%(QUALITY, min_length, overlap, infofile,\
                 FWDPRIMER, REVPRIMER, REV_COMPLEM_REVPRIMER, REV_COMPLEM_FWDPRIMER, outfile1, outfile2, read1, read2) \
                 #FAD_FP, RAD_RP, REV_COMPLEM_RAD_RP, REV_COMPLEM_FAD_FP, outfile1, outfile2, read1, read2)
    output = Popen(cmd, shell=True, stdout=PIPE).communicate()[0]
    #print cmd, '\n', output
    return output


'''
cutadapt -q 30 --minimum-length 100 --overlap 10 -g fwd_adapter_primer -G rev_adapter_primer -a rev_complem_rev_adapter_primer -A rev_complem_fwd_adapter_primer -n 3 -o out.1.fastq -p out.2.fastq reads.1.fastq reads.2.fastq

-g removes (i)
-G removes (iii)
-a removes (ii) , so give the rev complem.
-A removes (iv)
-n --times option. Cutadapt will then search for all the given adapter sequences repeatedly, either until no adapter match was found or until the specified number of rounds was reached. I think 2 is the max adapters I will see per read, but I still gave 3.

Although pandaseq also has the option to provide -p fwdprimerand -q revprimer, I am not sure if it is also removing (ii) and (iv). So for now I will do cutadapt and check adapter/primer removing ability of pandaseq later.

I will enforce the following Parameters in cutadapt:
- quality threshold (30)
- overlap of pattern (to remove) and read (10)
- keep both trimmed and untrimmed reads (trim but do not discard)
- min length (100)
- output
- info file
- paired output

'''
'''


######## I want to avoid cases where it is trimming seqs of length 3 or which can be random
(high expect.), so use the --overlap criteria.


No. of allowed errors:
0-9 bp: 0; 10-19 bp: 1; 20-29 bp: 2; 30-39 bp: 3; 40-49 bp: 4; 50-59 bp: 5; 60-65 bp: 6

Overview of removed sequences (5')
length  count   expect  max.err error counts
3       22      97.3    0       22
4       9       24.3    0       9

Overview of removed sequences (3' or within)
length  count   expect  max.err error counts
3       48      97.3    0       48
4       21      24.3    0       21
5       14      6.1     0       14

###### Also some seqs are probably in the middle and being treated as 3' causing almost the entire read
to be trimmed, so either use different parameter (-a or -g)

69      5       0.0     6       2 1 0 1 0 0 1
70      59      0.0     6       17 17 8 1 6 4 6
71      44      0.0     6       13 10 6 4 5 1 5
72      223     0.0     6       61 46 33 15 23 23 22
73      90      0.0     6       16 17 12 10 13 13 9
74      5       0.0     6       0 1 2 0 2
75      7       0.0     6       2 3 0 1 1
76      7       0.0     6       1 2 1 1 0 1 1
77      9       0.0     6       4 1 0 1 0 2 1
78      4       0.0     6       0 1 0 1 1 0 1
79      1       0.0     6       0 0 1
80      6       0.0     6       0 2 0 0 2 1 1
81      3       0.0     6       0 1 0 0 1 1
82      5       0.0     6       1 2 0 0 0 0 2
83      8       0.0     6       1 2 3 0 0 1 1
84      3       0.0     6       0 0 0 1 1 0 1
85      2       0.0     6       0 1 0 0 0 0 1
88      1       0.0     6       1
93      1       0.0     6       0 1
94      1       0.0     6       0 1
100     1       0.0     6       0 1
103     1       0.0     6       0 1
106     2       0.0     6       1 0 0 0 1
108     1       0.0     6       1
110     1       0.0     6       0 1
111     1       0.0     6       0 1
115     1       0.0     6       0 1
121     1       0.0     6       1
122     1       0.0     6       1
125     1       0.0     6       1
133     2       0.0     6       1 0 0 0 0 0 1
147     1       0.0     6       1
152     1       0.0     6       1
157     1       0.0     6       1
158     1       0.0     6       1
182     1       0.0     6       1
'''

if __name__=='__main__':
    datetime = strftime("%a, %d %b %Y %I:%M:%S %p", localtime())
    cmd = 'echo ' + datetime 
    os.system(cmd)
    main(sys.argv)

