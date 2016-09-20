import os
import sys
from utils import *
import operator
from time import localtime, strftime
import argparse
import re

# Usage: python ~/scripts/python/renumber_new_otus.py -i bisa-signif-taxa-pval-0.01-iv-70-remF8cd80k-table_mc80000.biom.txt -c /home/williamslab/qiime_software/python-2.7.3-release/lib/python2.7/site-packages/picrust/data/16S_13_5_precalculated.tab

def main(args):
    # calculate the following stats for bisa significant taxa:
    # number of otus under the 'level' (genus) (here, multiple occurances of the same species are treated as DIFFERENT otus)
    # number of distinct species (multiple occurances of the same species are treated as SAME otu)
    # number of seqs.
    parser = argparse.ArgumentParser(description='check which taxa names (if any) have different copy numbers using their original/standard otu ids.')
    parser.add_argument('-i', '--infile')
    parser.add_argument('-o', '--outfile', default="renumbered-non-new-otu-ids.txt")  # output filename
    parser.add_argument('-g', '--groupsfile')  # file containing the group info for each sample
    parser.add_argument('-c', '--categories' , nargs='*')  \
              # creates a list of groups (categories) from the group file to calc. stats for, currently only one category
    parser.add_argument('-d', '--delimiter', default='\t')  # delimiter for file
    parser.add_argument('-l', '--leveltosplit', default='s')  # taxonomic level to split. options allowed: p,c,o,f,g,s
    parser.add_argument('-n', '--otuidcolumn', default=0, type=int) # 0-index based column used for otu ids
    parser.add_argument('-t', '--taxonnamecolumn', default=-1, type=int) # 0-index based column used for taxonomy names
    
    args = parser.parse_args()

    if len(sys.argv)==1 :
        parser.print_help()
        sys.exit('\natleast one argument required\n')

    infile  = args.infile
    outfile = args.outfile
    delim = args.delimiter
    old_otu_ids_col = args.otuidcolumn
    taxon_col = args.taxonnamecolumn
    tax_level_to_split = args.leveltosplit
    
    categories_to_print = ['Group']
    if args.categories != None:
        categories_to_print = args.categories

    group_file = ''
    if args.groupsfile != None:
        group_file = args.groupsfile
    mapping_file = read_file(group_file) 

    for c in categories_to_print:
        check_infile(c, infile,outfile, delim,old_otu_ids_col,taxon_col, mapping_file, tax_level_to_split)

def check_infile(c, infile,outfile, delim,old_otu_ids_col,taxon_col, mapping_file,tax_level_to_split):
    # find index of SampleID and category to be summarized
    labels = mapping_file[0].split(delim)
    indx_sampleid = indx_categ = ''

    try:
        indx_sampleid = labels.index("#SampleID")
    except ValueError:
        print "SampleID not in list."

    try:
        indx_categ = labels.index(c)
    except ValueError:
        print "%s not in list." %c

    categ_samples_dict = list_to_dict(mapping_file, delim, ',', "current", indx_categ, indx_sampleid)
    category_list_sorted = sorted(categ_samples_dict) # so the categories are sorted (alphabetically)


    d_list=list()
    comment_list=list()
    # collect samples belonging to the same group
    categ_col_indx_dict = dict() # 0 based indexes, key is category, value is cols of samples in that category
 
    for l in open(infile):
        l = l.strip()
        if '#' not in l:
            d_list.append(l)
        else:
            comment_list.append(l)
            if '#OTU ID' in l or 'taxonomy' in l:
                array = l.split(delim)
                if old_otu_ids_col == array.index('#OTU ID'):
                    pass
                else:
                    print array.index('#OTU ID')
                    sys.exit('Default OTU id column number is incorrect. Use -n \n')
                if array[taxon_col] == 'taxonomy': # using index would fail for taxon, koz default is -1
                    pass
                else:
                    print array.index('taxonomy')
                    sys.exit('Default taxonomy column number is incorrect. Use -t \n')
                # now you know correct ids of otu id and taxonomy               
                for k in category_list_sorted:
                   v = categ_samples_dict[k]
                   indxs = [array.index(i) for i in v.split(',') if i in array]
                   categ_col_indx_dict[k] = indxs
                continue

    do_something(d_list, tax_level_to_split,category_list_sorted)

def do_something(c, mapping_file, infile, outpfile, delim):
    # number of otus under the 'level' (genus) (here, multiple occurances of the same species are treated as DIFFERENT otus)
    # number of distinct species (multiple occurances of the same species are treated as SAME otu)
    # number of seqs.
    same_taxon_same_otu = dict() # number of distinct species (multiple occurances of the same species are treated as SAME otu)
    same_taxon_diff_otus = dict() # list of otus under the 'level' (genus) (here, multiple occurances of the same species are treated as DIFFERENT otus)
    #same_taxon_diff_otus_std_err = dict()
    # I1  I2  I3  N1  N2  N3  taxon
    #  2   4  1    0   1  2   g1_s1
    #  1   3  3    1   0  3   g1_s1
    #  1   3  3    1   0  3   g1_s2
    ### for I1, genus g1 has 2 distinct species (s1,s2), 3 otus (s1,s1,s2)
    ### for I2, genus g1 has 2 distinct species (s1,s2), 3 otus (s1,s1,s2)
    ### for I3, genus g1 has 2 distinct species (s1,s2), 3 otus (s1,s1,s2)
    ### for I, genus g1 has 2 distinct species (s1,s2), 3 otus (avg. s1 + avg. s1 + avg. s2)
    same_taxon_occurance_counter = dict() #counter for occurances of otus with same taxon
    for l in d_list:
        contents = l.split(delim)
        # need a file of patterns to keep track of genus under consideration or split the line with '; s__' or ': s__'
        taxon = contents[taxon_col] # includes species level info.
        taxon_name = ''  # genus level name
        pattern_to_split1 = ';' + ' ' + tax_level_to_split + '__' # '; s__'
        pattern_to_split2 = ':' + ' ' + tax_level_to_split + '__' # ': s__'
        if pattern_to_split1 in taxon:
            taxon_name = taxon.split(pattern_to_split1)
        elif pattern_to_split2 in taxon:
            taxon_name = taxon.split(pattern_to_split2)
        temp_same_taxon_diff_otus = list()
        for k in category_list_sorted:
            v = categ_col_indx_dict[k]
            vals = [float(contents[i]) for i in v]
            #group_values_dict[k] = vals
            presence = []
            for i in vals:
                if i > 0: # only count the otu if there is atleast 1 seq in the sample
                    presence.append(1)
                else:
                    presence.append(0)
            #print taxon_name, k, v, vals
            temp_same_taxon_diff_otus[taxon_name].append(numpy.mean(presence))
        if taxon_name not in same_taxon_diff_otus:
            same_taxon_diff_otus[taxon_name] = [i for i in temp_same_taxon_diff_otus[taxon_name]]
        else:
            same_taxon_diff_otus[taxon_name] = [x+y for x,y in zip(same_taxon_diff_otus, temp_same_taxon_diff_otus)]
            #same_taxon_diff_otus[taxon_name] = [same_taxon_diff_otus[taxon_name][i] + temp_same_taxon_diff_otus[taxon_name][i] \
            #                      for i in range(len(temp_same_taxon_diff_otus[taxon_name]))]
               
################ to do from here ########################
    taxon_values_mean = dict() # 0 based indexes, key is taxon, value is tuple of means of values of samples in same category i.e. (mean_categ1, mean_categ2, mean_categ3)
    taxon_values_stderr = dict() # 0 based indexes, key is taxon, value is tuple of std err of values of samples in same category i.e. (err_categ1, err_categ2, err_categ3)



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


                    
'''
    taxon_otuid_dict = dict() # dict with taxon as key and (original) otu ids as value
    # print required info for topk
    for l in d_list:
         #print l
         # split contents of line
         array = l.split(delim)
         t = array[taxon_col]
         o = array[old_otu_ids_col]
         if t not in taxon_otuid_dict:
             taxon_otuid_dict[t] = list()
         taxon_otuid_dict[t].append(o)

    otuid_copy_numb_dict = list_to_dict(read_file(copynumb_file)) # dict with (original) otu ids as key and copy number as value

    taxon_copy_numb_dict = dict() # dict with taxon as key and copy numbers as value
    for t in taxon_otuid_dict:
        if t not in taxon_copy_numb_dict:
            taxon_copy_numb_dict[t] = list()
        otus = taxon_otuid_dict[t]
        for o in otus:
            if o in otuid_copy_numb_dict:
                c = otuid_copy_numb_dict[o]
                taxon_copy_numb_dict[t].append(c)
            else: # these are 'New*' otu ids, so absent in the gg 16s count files
                pass

    count = 0
    for t in taxon_otuid_dict:
        v = set(taxon_copy_numb_dict[t])
        if len(v) != 1:
            print t, v
            count += 1
    print '%d taxa, whose otu ids have different copy numbers' %count

'''
#    writeLIST_to_file(list, outfile)

  

if __name__=='__main__':
    datetime = strftime("%a, %d %b %Y %I:%M:%S %p", localtime())
    cmd = 'echo ' + datetime 
    os.system(cmd)
    main(sys.argv)

