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

#usage: python ~/scripts/python/randomize_core_microbiome.py -i /home/williamslab/Desktop/test/table_mc37500.biom -g /home/williamslab/Desktop/test/swg-sample-sheet-qiime-controls-v2-e3.txt -c Cultivar -n 5


def main(args):

    parser = argparse.ArgumentParser(description='Randomize data and calc. significance')
    parser.add_argument('-i', '--infile')
    parser.add_argument('-o', '--outfile', default="results")  # output filename
    parser.add_argument('-c', '--categories' , nargs='*') \
        # creates a list of categories to be summarized (enter as space delimited on cmd line), currently only one category
    parser.add_argument('-g', '--groupsfile')  # file containing the group where each sample belongs
    parser.add_argument('-f', '--pvaladjmethod', default='bf')  # pval adjust method: bonferroni (bf), benjamini hochberg (bh) 
    parser.add_argument('-d', '--delimiter', default='\t')  # delimiter for file
    parser.add_argument('-t', '--imagetype', default='pdf')  # file type of images, default pdf, other allowed are png or jpg
    parser.add_argument('-n', '--ntimes', default=1000, type=int) # randomize it ntimes

    args = parser.parse_args()

    mpl.rc('font', family="Arial", size=10) #http://matplotlib.org/users/customizing.html
    if len(sys.argv)==1 :
        parser.print_help()
        sys.exit('\natleast one argument required\n')

    infile  = args.infile
    outpfile = args.outfile
    delim = args.delimiter
    f_type = args.imagetype
    p_val_adj = args.pvaladjmethod
    categories_to_print = ['Group']
    if args.categories != None:
        categories_to_print = args.categories
    
    group_file = ''
    if args.groupsfile != None:
        group_file = args.groupsfile
    mapping_file = read_file(group_file)

    global logfile, logfilelist, NTIMES
    logfile = infile + '_' + '_'.join(categories_to_print) + '_' + 'core_microb-log-file.txt'
    logfilelist = list()
    NTIMES = args.ntimes

    infile_txt = ''
    # check if biom or txt
    if infile.endswith('.biom'):
        # convert biom to txt
        cmd = 'biom convert -i %s -o %s.txt -b --header-key taxonomy' %(infile, infile) # pattern
        output = Popen(cmd, shell=True, stdout=PIPE).communicate()[0]  # comment to avoid overqrite warning
        infile_txt = infile + '.txt'
    elif infile.endswith('.txt'):
        infile_txt = infile
    else:
        sys.exit('\ninfile can either be tab-delimited (.txt) or biom-format (.biom)\n')
    lines = read_file(infile_txt)
    if 'OTU ID' not in lines[1] and 'taxonomy' not in lines[1] : ### see what this should be ## 
        sys.exit('\nheader line is required, tab-delimited, contains taxonomy and samples\n')
    
    for c in categories_to_print:
        calc_significance(c, mapping_file, lines, infile, outpfile, delim, f_type, p_val_adj, group_file)

'''
* convert text to biom format
biom convert -i otu_table_mc2_w_tax_no_pynast_failures.biom.txt -o otu_table_mc2_w_tax_no_pynast_failures_biom_to_txt_to_biom.biom --table-type="otu table"
'''

'''
Randomize:
 the group (e.g. A, D) of the sample....same as randomizing the sample label (from the otu table), so it will get different column of values
 The otu label, so it will get different row of values
 The numbers in the otu table.
'''

def shuffle_dict(a):
    keys = a.keys()
    values = list()
    lengths = list()
    for i in a.values():
        v = i.split(',')
        values.extend(v)
        lengths.append(len(v)) # sizes of original lists
    random.shuffle(values)
    new_values = divide_list(values, lengths) # divide the shuffled values as per original sizes
    #print values, '\n', new_values
    return dict(zip(keys, new_values))



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



def compile_results(core_otus_file, delim): # get unique elements from last column (otus) from the core_otu_100 files
    taxon_list = list()
    result_ = read_file(core_otus_file)
    for l in result_:
        l = l.strip()
        contents = l.split(delim)
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



def calc_significance(c, mapping_file, lines, infile, outpfile, delim, f_type, p_val_adj, group_file):
    # find index of SampleID and category to be summarized  # e.g. swg or non-swg
    labels = mapping_file[0].split(delim)
    indx_sampleid = indx_categ = ''

    try: # this is from the mapping file
        indx_sampleid = labels.index("#SampleID")
    except ValueError:
        print "SampleID not in list."

    try:
        indx_categ = labels.index(c)
    except ValueError:
        print "%s not in list." %c

    group = 'A'  ### to do : change this to 'swg'
    print 'Running true data'
    exec_core_microb_cmd(infile, 'true_result', group_file, c, group) # run core microbiome cmd on original data

    frac_samples = 100
    true_result_file_frac_samples = './%s/core_otus_%s.txt' %('true_result', frac_samples)
    true_results_otu_list = compile_results(true_result_file_frac_samples, delim)

    categ_samples_dict = list_to_dict(mapping_file, delim, ',', "current", indx_categ, indx_sampleid)
    #print categ_samples_dict

    # randomize n times
    for i in range(NTIMES):
        print 'Running randomized data: %s' %i
        d = categ_samples_dict
        # so the categories are shuffled, i.e. randomly shuffle what keys correspond to what values
        categ_samples_dict_shuffled = shuffle_dict(d)
        #print categ_samples_dict_shuffled

        # write the shuffled group to file
        new_mapping_file = outpfile+'-mapping-file.txt'
        write_shuffled_dict_new_file(categ_samples_dict_shuffled, c, new_mapping_file)

        exec_core_microb_cmd(infile, i, new_mapping_file, c, group) # run core microbiome cmd.

    # compile the results
    all_results_otu_list = list()
    for i in range(NTIMES):
        result_file_frac_samples = './%s/core_otus_%s.txt' %(i, frac_samples)
        taxons_ = compile_results(result_file_frac_samples, delim) # unique set of core taxa from randomized data
        all_results_otu_list.extend(taxons_)

    # calculate stats after completing the for loop.
    randomized_otus = calc_freq_elem_list(all_results_otu_list) # dict of otus and freq of occurance from random data
    #print randomized_otus 
    writeDICT_to_file(randomized_otus, outpfile+'-randomized-mapping-file-results-'+str(frac_samples)+'.txt')

    # print the number of random occurances for each true core microbiome otu
    for o in true_results_otu_list:
        if o in randomized_otus:
            freq = int(randomized_otus[o])
            if freq/float(NTIMES) < 0.05:
                print o , freq

    ### some cmds to delete the random data folders ### 
    for i in range(NTIMES):
        cmd = 'rm -r %s' %(i) 
        output = Popen(cmd, shell=True, stdout=PIPE).communicate()[0]

    sys.exit(0)

'''
    for i in range(NTIMES):
        # shuffle otus
        otu_names_list_tup = readColumnsSep(f, '\t', 0)
        otu_names = [i[0] for i in otu_names_list_tup] # list of otu names
        counter = 0
        for otu in random.shuffle(otu_names):
            line = lines[counter]
            cols = line.split(delim)
            print otu , '\t', cols[1..n]  # fix this to do
            counter += 1
        # write the shuffled to file
        # run core microbiome cmd.
        # compile the results
    # calculate stats after completing the for loop.
'''
    
'''    
    # collect samples belonging to the same group
    categ_col_indx_dict = dict() # 0 based indexes, key is category, value is cols of samples in that category
    taxon_values_mean = dict() # 0 based indexes, key is taxon, value is tuple of means of values of samples in same category i.e. (mean_categ1, mean_categ2, mean_categ3)
    taxon_values_stderr = dict() # 0 based indexes, key is taxon, value is tuple of std err of values of samples in same category i.e. (err_categ1, err_categ2, err_categ3)
    for l in lines:
        #print l, "------------------"
        l = l.strip()
        contents = l.split(delim)
        if '#' in l or not l:
            continue
        if 'Taxon' in l:
            for k in category_list_sorted:
               v = categ_samples_dict[k]
               indxs = [contents.index(i) for i in v.split(',') if i in contents]
               categ_col_indx_dict[k] = indxs
            continue
        taxon_name = contents[0]
        taxon_values_mean[taxon_name] = list()
        taxon_values_stderr[taxon_name] = list()
        group_values_dict = dict()
        # calc. mean, std. err for each group
        for k in category_list_sorted:
            v = categ_col_indx_dict[k]
            vals = [float(contents[i]) for i in v]
            group_values_dict[k] = vals
            #print taxon_name, k, v, vals
            taxon_values_mean[taxon_name].append(str(numpy.mean(vals)))  # convert these directly to float
            taxon_values_stderr[taxon_name].append(str(numpy.std(vals)/math.sqrt(len(vals))))
        #print taxon_values_mean, taxon_values_stderr
        two_sample_stats_test(taxon_name, group_values_dict, contrasts, p_val_adj, variance, delim, len(lines)-1, stats_test) # the first line is header
    #print taxon_values_mean, taxon_values_stderr
    draw_stacked_bar_w_err_bar(c, taxon_values_mean, taxon_values_stderr, category_list_sorted, outpfile, f_type)
    writeLIST_to_file(logfilelist, logfile)
'''



if __name__=='__main__':
    datetime = strftime("%a, %d %b %Y %I:%M:%S %p", localtime())
    cmd = 'echo ' + datetime 
    os.system(cmd)
    main(sys.argv)

