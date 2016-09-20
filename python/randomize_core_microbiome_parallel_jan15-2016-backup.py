import os
import sys
sys.path.insert(0,"/usr/lib/python2.7/dist-packages/") 
from utils import *
import operator
from time import localtime, strftime
import argparse
import math
import numpy
import matplotlib.pyplot as plt
import brewer2mpl
import pylab
from scipy import stats
#from statsmodels.stats import multicomp
import matplotlib as mpl
import string
import random
from subprocess import Popen, PIPE
import collections
from multiprocessing import Pool, Lock
import itertools

#usage: python randomize_core_microbiome_parallel.py -i table_mc37500.biom -m swg-sample-sheet-qiime-controls-v2-e3.txt -c Cultivar -n 4 -g A -f 100



def main(args):

    parser = argparse.ArgumentParser(description='Randomize data and calc. significance of core microbiome')
    parser.add_argument('-i', '--infile')
    parser.add_argument('-o', '--outfile', default="results")  # output filename
    parser.add_argument('-c', '--categories' , nargs='*') \
        # creates a list of categories to be summarized (enter as space delimited on cmd line), currently only one category
    parser.add_argument('-m', '--mappingfile')  # file containing the group where each sample belongs i.e. maps each sample to a group
    parser.add_argument('-g', '--group', default='swg')  # group e.g A,D,swg for whom the core microbiome is to be calculated
    parser.add_argument('-p', '--pvaladjmethod', default='bf')  # pval adjust method: bonferroni (bf), benjamini hochberg (bh) 
    parser.add_argument('-d', '--delimiter', default='\t')  # delimiter for file
    parser.add_argument('-t', '--imagetype', default='pdf')  # file type of images, default pdf, other allowed are png or jpg
    parser.add_argument('-n', '--ntimes', default=1000, type=int) # randomize it ntimes
    parser.add_argument('-f', '--fracsamples', default=100, type=int) # otu present in atleast fraction of samples (give values 50-100)

    args = parser.parse_args()

    #mpl.rc('font', family="Arial", size=10) #http://matplotlib.org/users/customizing.html
    if len(sys.argv)==1 :
        parser.print_help()
        sys.exit('\natleast one argument required\n')

    global FRAC_SAMPLES, DELIM, GROUP, NTIMES, OUTPFILE

    infile  = args.infile

    NTIMES = args.ntimes
    OUTPFILE = args.outfile
    DELIM = args.delimiter
    FRAC_SAMPLES = args.fracsamples

    f_type = args.imagetype
    p_val_adj = args.pvaladjmethod

    categories_to_print = ['Group']
    if args.categories != None:
        categories_to_print = args.categories
    
    group_file = ''  # name of mapping (group) file
    if args.mappingfile != None:
        group_file = args.mappingfile
    mapping_info_list = read_file(group_file) # reads the mapping (group) info as a list

    GROUP = args.group

    global logfile, logfilelist, INFILE
    logfile = infile + '_' + '_'.join(categories_to_print) + '_' + 'core_microb-log-file.txt'
    logfilelist = list()
    INFILE = infile

    infile_txt = ''
    # check if biom or txt
    if infile.endswith('.biom'):
        # convert biom to txt
        cmd = 'biom convert -i %s -o %s.txt -b --header-key taxonomy' %(infile, infile)
        # newer biom version need something like: biom convert -i table_mc37500.biom -o table_mc37500.biom.txt --header-key taxonomy --to-tsv
        output = Popen(cmd, shell=True, stdout=PIPE).communicate()[0]
        infile_txt = infile + '.txt'
    elif infile.endswith('.txt'):
        infile_txt = infile
        # convert text to biom format
        #biom convert -i otu_table_mc2_w_tax_no_pynast_failures.biom.txt -o otu_table_mc2_w_tax_no_pynast_failures_biom_to_txt_to_biom.biom --table-type="otu table"
    else:
        sys.exit('\ninfile can either be tab-delimited (.txt) or biom-format (.biom)\n')

    lines = read_file(infile_txt)
    if 'OTU ID' not in lines[1] and 'taxonomy' not in lines[1] :
        sys.exit('\nheader line is required, tab-delimited, contains taxonomy and samples\n')
    
    for c in categories_to_print:
        calc_significance(c, mapping_info_list, lines, f_type, p_val_adj, group_file)

'''
Randomize:
 the group (e.g. A, D) of the sample....same as randomizing the sample label (from the otu table), so it will get different column of values
 The otu label, so it will get different row of values
 The numbers in the otu table.
'''

def shuffle_dict(iternumb): ## takes iternumb which is number of random iteration
    a = CATEG_SAMPLES_DICT
    keys = a.keys()
    values = list()
    lengths = list()
    for i in a.values():
        v = i.split(',')
        values.extend(v)
        lengths.append(len(v)) # sizes of original lists
    random.shuffle(values)
    new_values = divide_list(values, lengths) # divide the shuffled values as per original sizes
    categ_samples_dict_shuffled = dict(zip(keys, new_values))

    # write the shuffled group to file
    new_mapping_file = OUTPFILE + '-' + str(iternumb) + '-mapping-file.txt'  #### renamed this so it writes to a new file and we can parallelize
    write_shuffled_dict_new_file(categ_samples_dict_shuffled, CATEGO, new_mapping_file) 
    exec_core_microb_cmd(INFILE, iternumb, new_mapping_file, CATEGO, GROUP) # run core microbiome cmd.
    return



def divide_list(a, lengths):
    return a[:lengths[0]] , a[lengths[0]:]  # a[start:end] # items start through end-1



def write_shuffled_dict_new_file(DICT, categ, filename):
    outfile = open(filename,'w') # make this 'a' to append while debugging
    outfile.write("%s\t%s\n" % ('#SampleID' , categ))
    for k, v in DICT.items():
        for i in v:
            outfile.write("%s\t%s\n" % (i , k))
    outfile.close()
    print 'DICT written to file %s' % (filename)
    return



def compile_results(iternumb): # get unique elements from last column (otus) from the core_otu_100 files
    core_otus_file = './%s/core_otus_%s.txt' %(iternumb, FRAC_SAMPLES)
    taxon_list = list()
    result_ = read_file(core_otus_file)
    for l in result_:
        l = l.strip()
        contents = l.split(DELIM)
        if '#' in l or not l:
            continue
        taxon_list.append(contents[1])
    return list(set(taxon_list))



def calc_freq_elem_list(a):
    counter=collections.Counter(a)
    return counter



def exec_core_microb_cmd(infile, i, new_mapping_file, categ, group):
    cmd = 'compute_core_microbiome.py -i %s -o %s --mapping_fp %s --valid_states "%s:%s"' %(infile, i, new_mapping_file, categ, group)
    output = Popen(cmd, shell=True, stdout=PIPE).communicate()[0]    
    return



def calc_significance(c, mapping_info_list, lines, f_type, p_val_adj, group_file):
    # find index of SampleID and category to be summarized  # e.g. swg or non-swg
    labels = mapping_info_list[0].split(DELIM)
    indx_sampleid = indx_categ = ''

    try: # this is from the mapping file
        indx_sampleid = labels.index("#SampleID")
    except ValueError:
        print "SampleID not in list."

    try:
        indx_categ = labels.index(c)
    except ValueError:
        print "%s not in list." %c

    global CATEGO
    CATEGO = c

    # run core microbiome on original data and compile results
    print 'Running true data'
    original_data = 'true_result'
    exec_core_microb_cmd(INFILE, original_data, group_file, CATEGO, GROUP)
    true_results_otu_list = compile_results(original_data)

    categ_samples_dict = list_to_dict(mapping_info_list, DELIM, ',', "current", indx_categ, indx_sampleid)

    global CATEG_SAMPLES_DICT
    CATEG_SAMPLES_DICT = categ_samples_dict
    
    pool = Pool()              # 1 process per core
    #pool = Pool(processes=2)              # process per core

    # randomize n times
    print 'Running randomized data:'
    pool.map(shuffle_dict, range(NTIMES))  # do shuffle_dict n times with pool
    
    # compile the results from randomization
    taxons_ = pool.map(compile_results, range(NTIMES))  # this returns a list of list i.e. collects the unique set of core taxa from each randomized data
    all_results_otu_list = [item for sublist in taxons_ for item in sublist] #  taxons_ -> sublist -> item

    # calculate freq of otu being a core microb from the randomizations
    randomized_otus = calc_freq_elem_list(all_results_otu_list) # dict of otus and freq of occurance from random data
    writeDICT_to_file(randomized_otus, OUTPFILE+'-randomized-mapping-file-results-'+str(FRAC_SAMPLES)+'.txt')

    # print the number of random occurances for each true core microbiome otu (checks significance)
    signif_core_microb_otu_list = list()
    for o in true_results_otu_list:
        if o in randomized_otus:
            freq = int(randomized_otus[o])
            if freq/float(NTIMES) < 0.05:
                otu = '%s\t%s' % (o, freq)
                print otu
                signif_core_microb_otu_list.append(otu)

    writeLIST_to_file(signif_core_microb_otu_list, OUTPFILE+'-signif_core_microb_otu_list-results-'+str(FRAC_SAMPLES)+'.txt')


    ### some cmds to delete the randomly generated data folders and mapping files ### 
    for i in range(NTIMES):
        cmd = 'rm -r %s' %(i) 
        output = Popen(cmd, shell=True, stdout=PIPE).communicate()[0]
        cmd = 'rm ' + OUTPFILE + '-' + str(i) + '-mapping-file.txt'
        output = Popen(cmd, shell=True, stdout=PIPE).communicate()[0]
    sys.exit(0)



if __name__=='__main__':
    datetime = strftime("%a, %d %b %Y %I:%M:%S %p", localtime())
    cmd = 'echo ' + datetime 
    os.system(cmd)
    main(sys.argv)

