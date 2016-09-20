import os
import sys
from utils import *
import operator
from time import localtime, strftime
import argparse
import re
import os.path
from subprocess import Popen, PIPE


#usage: 
# the following was run in Trimmed_fastq_files/q30
# python ~/scripts/python/stitch_pe_reads_pandaseq_nextera.py -i ../../../../analysis/michigan-illumina.txt -q -o ../../pandaseq_stitched/q30


def main(args):

    parser = argparse.ArgumentParser(description='Check if the barcodes  are present in the sample')
    parser.add_argument('-i', '--infile') # file containing the sample name and the barcode string
    parser.add_argument('-o', '--outdir', default="./out")  # output dirname
    parser.add_argument('-l', '--logdir', default="./log")  # log dirname
    parser.add_argument('-u', '--unpaireddir', default="./unpaired_files/")  # unpaired files dirname
    parser.add_argument('-f', '--fasta', action='store_true', default=False) # files in fasta format, default fastq
    parser.add_argument('-q', '--qiimeformat', action='store_true', default=False) # output files in qiime compatible fasta format, default pandaseq stitched fastq format
    parser.add_argument('-d', '--delimiter', default='\t')  # delimiter for file

    args = parser.parse_args()

    infile  = args.infile
    outdir = args.outdir if args.outdir[-1] == '/' else args.outdir + '/'
    logdir = args.logdir if args.logdir[-1] == '/' else args.logdir + '/'
    unpaired_dir = args.unpaireddir if args.unpaireddir[-1] == '/' else args.unpaireddir + '/'
    filetype = '.fasta' if args.fasta else '.fastq'
    qiimeformat = args.qiimeformat
    delim = args.delimiter

    SAMPLE_BARCODE_DICT = dict()
    sampleid_col = fbarcode_col = rbarcode_col = 0

    if qiimeformat: 
        if infile == None :
            parser.print_help()
            sys.exit('\nTo output files in qiime format, file with sample and barcode is required\n')
        else:
            lines = read_file(infile)
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
                barcode = cont[fbarcode_col] + cont[rbarcode_col] # a '_' separating these throws an error in mapping file validation
                if sample in SAMPLE_BARCODE_DICT:
                    sys.exit('\nduplicate samples/barcodes detected, check sample sheet\n')
                else:
                    SAMPLE_BARCODE_DICT[sample] = barcode

    for dir in [outdir, logdir, unpaired_dir]:
        if not os.path.exists(dir):
            os.makedirs(dir)

    global R_1, R_2
    R_1 = '_L001_R1_001'
    R_2 = '_L001_R2_001'
    R_12 = '_L001_R12_001'
    
    completed_files_list = list()
    for f in sorted(os.listdir('./')):
        if os.path.isfile(f) and f.endswith(filetype):
            for k in SAMPLE_BARCODE_DICT:  # for each element in dict
                if f.startswith('trim_'+k) and f not in completed_files_list:   # if the key is prefix in filename
                    #put the read 1 and read 2 file into completed list
                    #go ahead and do the stitching
                    barcode = SAMPLE_BARCODE_DICT[k]
                    read1 = read2 = ''
                    if R_1 in f:
                        read1 = f
                        read2 = f.replace(R_1, R_2)
                    elif R_2 in f:
                        read2 = f
                        read1 = f.replace(R_2, R_1)
                    completed_files_list.extend([read1, read2])
                    tmp = read1.replace(R_1, '_log')
                    logfile = logdir + tmp.replace(filetype, '.txt')
                    tmp = outdir + read1.replace(R_1, R_12)
                    outputfile = cmd = ''
                    tmp1 = read1.replace(R_1, '_unpaired')
                    unpaired_file = unpaired_dir + tmp1.replace(filetype, '.txt')
                    if filetype == '.fasta':
                        outputfile = tmp.replace(filetype, '.fasta')
                        cmd='pandaseq -f %s -r %s -B -g %s -U %s -w %s' %(read1, read2, logfile, unpaired_file, outputfile)
                    else:
                        outputfile = tmp
                        cmd='pandaseq -f %s -r %s -B -F -g %s -U %s -w %s' %(read1, read2, logfile, unpaired_file, outputfile)
                    #print cmd
                    #sys.exit(0)
                    output_start = Popen(cmd, shell=True, stdout=PIPE).communicate()[0]
                    if qiimeformat and filetype != '.fasta':
                        sname = k+'_R12_'
                        outsname = outputfile
                        outputfile = outdir + sname
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

