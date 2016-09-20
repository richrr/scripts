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



#usage: python ~/scripts/python/check_barcodes.py -i ../../../../../SampleSheet_w_groups.csv -d , -o barcode_stats.csv

#python ~/scripts/python/check_barcodes.py -i ../../../../../SampleSheet_w_groups.csv -d , -o barcode_stats_q25.csv -l its-log-q25.txt -q 25 -x ./Intermed_fastq_files/q25/ -y ./Intermed_info_files/q25/ -z ./Trimmed_fastq_files/q25/

# for using only read 2, uncomment (A-E) and comment lines directly above them if present (except for D and E)
#python ~/scripts/python/check_barcodes.py -i ../../../../../SampleSheet_w_groups.csv -d , -o barcode_stats_read2_q30.csv -l its-log-read2_q30.txt -x ./Intermed_fastq_files/read2_q30/ -y ./Intermed_info_files/read2_q30/ -z ./Trimmed_fastq_files/read2_q30/
    

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

def grep_two_patterns(start, end, f):
    cmd = 'egrep -c "^%s.+%s" %s' %(start, end,f) # pattern
    output = Popen(cmd, shell=True, stdout=PIPE).communicate()[0]  # grep -c 'pattern' filename
    #print cmd , '\n', output
    '''
    cmd = 'egrep -c "^%s.+%s$" %s' %(start, end,f) # pattern
    output = Popen(cmd, shell=True, stdout=PIPE).communicate()[0]  # grep -c 'pattern' filename
    print cmd , '\n', output
    '''
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

def main(args):

    parser = argparse.ArgumentParser(description='Check if the barcodes  are present in the sample')
    parser.add_argument('-i', '--infile') # file containing the sample name and the barcode string
    parser.add_argument('-o', '--outfile', default="out.txt")  # output filename
    parser.add_argument('-l', '--logfile', default="ITS-log-file.txt")  # output filename
    parser.add_argument('-p', '--pairedend', action='store_true', default=False) # paired-end data, default single end
    parser.add_argument('-f', '--fasta', action='store_true', default=False) # files in fasta format, default fastq
    parser.add_argument('-q', '--quality', default=30, type=int) # Phred quality threshold of seq.
    parser.add_argument('-x', '--intermedfastq', default='./Intermed_fastq_files/')  # dir name for intermediate fastq files
    parser.add_argument('-y', '--intermedinfo', default='./Intermed_info_files/')  # dir name for intermediate info files
    parser.add_argument('-z', '--trimfastq', default='./Trimmed_fastq_files/')  # prefix for trimmed files
    #parser.add_argument('-g', '--groupsfile')  # file containing the group where each sample belongs
    #parser.add_argument('-c', '--insertcolumn', default=99999, type=int) # 0-index based column where the row_sums is to be printed
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
    #end_col = args.endcolum
    #insertcolumn = args.insertcolumn
    delim = args.delimiter
    if args.histogram:
        draw_histogram_adpaters_trimed_seqs(infile)
        sys.exit(0)

    global FWDPRIMER, FADAPTER, REVPRIMER, RADAPTER, QUALITY, R_1, R_2, INTERMFASTQ, INTERMINFO, TRIMFASTQF
    FWDPRIMER='ATCTACACTATGGTAATTGTGAACCWGCGGARGGATCA'
    FADAPTER='AATGATACGGCGACCACCGAG'       
    REVPRIMER='GCATCGATGAAGAACGCAGCGGCTGACTGACT'
    RADAPTER='ATCTCGTATGCCGTCTTCTGC'
    QUALITY=args.quality
    R_1 = '_L001_R1_001'
    R_2 = '_L001_R2_001'
    INTERMFASTQ = args.intermedfastq
    INTERMINFO = args.intermedinfo
    trimfastqd = args.trimfastq  # dir
    TRIMFASTQF = trimfastqd + 'trim_' # dir + prefix for files
    
    for dir in [INTERMFASTQ, INTERMINFO, trimfastqd]:
        if not os.path.exists(dir):
            os.makedirs(dir)
    
    lines = read_file(infile)
    #GRL4463_S49_L001_R1_001
    SAMPLE_BARCODE_DICT = dict()
    for l in lines:
        if '#' in l:
            continue
        cont = l.strip().split(delim)
        sample = cont[0] + '_S'
        barcode = cont[1]
        if cont[0] in SAMPLE_BARCODE_DICT:
            sys.exit('\nduplicate samples/barcodes detected, check sample sheet\n')
        else:
            SAMPLE_BARCODE_DICT[sample] = barcode
        # 5' ITS1FI2  forward primer 3'
        #       adapter                    fwd-primer
        #AATGATACGGCGACCACCGAG ATCTACACTATGGTAATTGTGAACCWGCGGARGGATCA
        #5'     adapter              barcode           rev-primer     3'
        #CAAGCAGAAGACGGCATACGAGAT GATTCCGGCTCA AGTCAGTCAGCCGCTGCGTTCTTCATCGATGC
        # reverse-complement
        #5'      rev-primer                barcode            adapter       3'
        #GCATCGATGAAGAACGCAGCGGCTGACTGACT TGAGCCGGAATC ATCTCGTATGCCGTCTTCTGCTTG
        # TGAGCCGGAATC
    outfile = open(outpfile, 'w')
    headers_list = ['file' , 'starts_bc', 'ends_bc', 'starts_rp', 'ends_rp', 'starts_ad', 'ends_ad', \
            'starts_bc_rad' , 'ends_bc_rad' ,  'starts_bc_ends_bc_rad', 'starts_bc_rad_ends_bc_rad',\
            'starts_rp_bc_rad' , 'ends_rp_bc_rad' ,  'starts_bc_ends_rp_bc_rad', 'starts_bc_rad_ends_rp_bc_rad',\
            'starts_rp_bc_rad_ends_rp_bc_rad',  'neg_starts_bc_ends_rp_bc_rad', 'neg_starts_bc_rad_ends_rp_bc_rad',\
             'starts_bc_neg_ends_rp_bc_rad', 'starts_bc_rad_neg_ends_rp_bc_rad', 'neg_starts_bc_rad_neg_ends_rp_bc_rad',\
             'starts_fp', 'ends_fp', 'contains_fp','starts_fad', 'ends_fad', 'contains_fad',\
             'starts_fad_fp', 'ends_fad_fp','contains_fad_fp']
    outfile.write ("%s\n" % delim.join(headers_list))                
    for f in sorted(os.listdir('./')):
        if os.path.isfile(f) and f.endswith(filetype):
            for k in SAMPLE_BARCODE_DICT:  # for each element in dict
                if k in f:   # if the key is prefix in filename
                #if k in f and R_2 in f:   # if the key is prefix in filename and only for read 2 --------------- A for read 2
                    barcode = SAMPLE_BARCODE_DICT[k]
                    counts_list = [f]
                    counts_list.extend(grep_string(barcode, f))
                    counts_list.extend(grep_string(REVPRIMER, f))
                    counts_list.extend(grep_string(RADAPTER, f))
                    bc_rad = barcode+RADAPTER
                    counts_list.extend(grep_string(bc_rad, f))
                    counts_list.append(grep_two_patterns(barcode, bc_rad, f))
                    counts_list.append(grep_two_patterns(bc_rad, bc_rad, f))
                    rp_bc_rad = REVPRIMER+barcode+RADAPTER
                    counts_list.extend(grep_string(rp_bc_rad, f))
                    counts_list.append(grep_two_patterns(barcode, rp_bc_rad, f))
                    counts_list.append(grep_two_patterns(bc_rad, rp_bc_rad, f))
                    counts_list.append(grep_two_patterns(rp_bc_rad, rp_bc_rad, f))
                    counts_list.append(grep_two_patterns_negate_first(barcode, rp_bc_rad, f))
                    counts_list.append(grep_two_patterns_negate_first(bc_rad, rp_bc_rad, f))
                    counts_list.append(grep_two_patterns_negate_second(barcode, rp_bc_rad, f))
                    counts_list.append(grep_two_patterns_negate_second(bc_rad, rp_bc_rad, f))
                    counts_list.append(grep_two_patterns_negate_both(bc_rad, rp_bc_rad, f))
                    
                    counts_list.extend(grep_string(FWDPRIMER, f))
                    counts_list.append(grep_string_contains(FWDPRIMER, f))
                    counts_list.extend(grep_string(FADAPTER, f))
                    counts_list.append(grep_string_contains(FADAPTER, f))
                    fad_fp = FADAPTER+barcode
                    counts_list.extend(grep_string(fad_fp, f))
                    counts_list.append(grep_string_contains(fad_fp, f))
                    str_counts_list = [str(i.strip()) for i in counts_list]
                    #print '\t'.join(str_counts_list)
                    outfile.write ("%s\n" % delim.join(str_counts_list))
                    #sys.exit(0)
    outfile.close()
    
    completed_files_list = list() 
    outfile = open(logfile, 'w')
    for f in sorted(os.listdir('./')):
        if os.path.isfile(f) and f.endswith(filetype):
            for k in SAMPLE_BARCODE_DICT:  # for each element in dict
                if k in f and f not in completed_files_list:   # if the key is prefix in filename
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


def remove_adap_bc_primer_quality_trim_min_length (read1, read2, barcode):
    # remove bc_rad from the front using -g and ^, 
    # rp_bc_rad from back using -b to account for partial seqs like pterSEQUENCE(of read1)
    global FWDPRIMER, FADAPTER, REVPRIMER, RADAPTER, QUALITY, R_1, R_2
    bc_rad    = barcode + RADAPTER
    rp_bc_rad = REVPRIMER + barcode + RADAPTER
    fad_fp    = FADAPTER + FWDPRIMER
    
    
    # remove adapter then do Quality trimming (and min length) separately
    # doing the above together causes quality trimm to be first
    # if the adapter, bc, rp bases are low quality, they might get trimmed
    # and then the trimming of adapter,bc,rp may be unsuccessful/confusing.
    # although in this case, following is needed only for read1, 
    # I am going to do to both just to be safe and for use of future datasets
    # although not needed, i am using the -a name=sequence format so i can use
    # the name for the adapter column in the info file
    output_dump = list()
    read_file = read_file_base = ''
    for read in [read1,read2]:
    #for read in [read2]:       #--------------- B for read 2
        output_dump.append('--- Read ---')
        adap_name = 'bc_rad__'
        param = '-g %s=^%s' %(adap_name, bc_rad) # --front, ^ enforces search at start of string
        output, infile_base, infile = trim_adap(adap_name, bc_rad, param, read, read)
        #cmd = 'cutadapt -g %s=^%s -o %s -e 0.05 -O %s --info-file=%s %s' %(adap_name, bc_rad, tmp, overlap, infofile, read)
        output_dump.append(output)
        
        adap_name = 'rp_bc_rad__'
        param = '-b %s=%s' %(adap_name, rp_bc_rad) # --anywhere
        output, infile_base, infile = trim_adap(adap_name, rp_bc_rad, param, infile_base, infile)
        #cmd = 'cutadapt -b %s=%s -o %s -e 0.05 -O %s --info-file=%s %s' %(adap_name, rp_bc_rad, tmp2, overlap, infofile, tmp)
        output_dump.append(output)    
        
        adap_name = 'fad_fp__'
        param = '-b %s=%s' %(adap_name, fad_fp) # --anywhere
        #param = '-q %d --minimum-length 100 -b %s=%s' %(QUALITY, adap_name, fad_fp)  #--------------- C for read 2
        output, infile_base, infile = trim_adap(adap_name, fad_fp, param, infile_base, infile)
        #cmd = 'cutadapt -b %s=%s -o %s -e 0.05 -O %s --info-file=%s %s' %(adap_name, fad_fp, tmp3, overlap, infofile, tmp2)
        output_dump.append(output)
        read_file = infile
        read_file_base = infile_base
    
    output1, output2 = min_length_filter_paired(read_file, read_file_base, QUALITY, R_1, R_2)         #--------------- D for read 2
    output_dump.extend([output1, output2])                                                            #--------------- E for read 2
    return output_dump 
       
def trim_adap(adap_name, adap, param, infile_base, infile):
    global INTERMFASTQ, INTERMINFO
    ofilename_base = adap_name + infile_base
    ofilename = INTERMFASTQ + ofilename_base
    infofile = INTERMINFO + ofilename_base + '_info.txt'
    overlap = len(adap) - 10
    cmd = 'cutadapt %s -o %s -e 0.05 -O %s --info-file=%s %s' %(param, ofilename, overlap, infofile, infile)
    output = Popen(cmd, shell=True, stdout=PIPE).communicate()[0]
    #print cmd, '\n', output
    return output, ofilename_base, ofilename

def min_length_filter_paired(read_file, read_file_base, QUALITY, R_1, R_2):
    global TRIMFASTQF
    trimmed2 = TRIMFASTQF + read_file_base  # read file 2
    trimmed1 = trimmed2.replace(R_2, R_1)
    read2file = read_file
    read1file = read_file.replace(R_2, R_1)
    cmd='cutadapt -q %d --minimum-length 100 --paired-output tmp.2.fastq -o tmp.1.fastq %s %s' %(QUALITY, read1file, read2file)
    output1 = Popen(cmd, shell=True, stdout=PIPE).communicate()[0]
    #print cmd, '\n', output1
    cmd='cutadapt -q %d --minimum-length 100 --paired-output %s -o %s tmp.2.fastq tmp.1.fastq' %(QUALITY, trimmed1, trimmed2)
    output2 = Popen(cmd, shell=True, stdout=PIPE).communicate()[0]
    #print cmd, '\n', output2
    cmd='rm tmp.1.fastq tmp.2.fastq'
    output = Popen(cmd, shell=True, stdout=PIPE).communicate()[0]
    return output1, output2

'''
#First trim the forward read, writing output to temporary files (we also add some quality trimming):
       cutadapt -q 10 --minimum-length 20 --paired-output tmp.2.fastq -o tmp.1.fastq reads.1.fastq reads.2.fastq
#Then trim the reverse read, using the temporary files as input:
       cutadapt -q 15 --minimum-length 20 --paired-output trimmed.1.fastq -o trimmed.2.fastq tmp.2.fastq tmp.1.fastq
#Remove the temporary files:
       rm tmp.1.fastq tmp.2.fastq
#quality threshold (30)
#overlap (length of pattern - 10)
#keep both trimmed and untrimmed reads (trim but do not discard)
#min length  (100)
#output
#info file
#paired output


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

