import os
import sys
from utils import *
import operator
from time import localtime, strftime
import argparse
import re
import numpy as np
import matplotlib
matplotlib.use('Agg') # to avoid server trying to display the figures
import matplotlib.pyplot as plt
import math

#[rodrrich@transkingdom Abx]$ python /nfs1/Morgun_Lab/richrr/scripts/python/intersection_plotter_overlap_selector.py -i /nfs1/Morgun_Lab/richrr/Abx/MiSeq1/analysis_read1/OTU_ClosedRef/miseq10_otu_table_filt_samp.biom.tsv.sorted.biom.CSS.tsv -j /nfs1/Morgun_Lab/richrr/Abx/MiSeq2/analysis_read1/OTU_ClosedRef/miseq11_otu_table.biom.tsv.sorted.biom.CSS.tsv -m /nfs1/Morgun_Lab/richrr/Abx/MiSeq1/miseq10_mapping.txt -n /nfs1/Morgun_Lab/richrr/Abx/MiSeq2/mapping.txt -g Treatment Experiment



def plot_JI_of_files(ListA, ListB, plot_outfile):
    #JIvalues is a list of pairs [x,y], where 'y' is the JI at cutoff 'x'
    JIvalues = jaccardIndex(ListA, ListB)
    plt.ioff()

    plt.figure(1)
    plt.subplot(211)
    plt.scatter(*zip(*JIvalues))
    plt.ylabel('Jaccard Index')

    JIvalues = jaccardIndex_counts(ListA, ListB)

    plt.subplot(212)
    plt.scatter(*zip(*JIvalues))
    plt.xlabel('Steps')
    plt.ylabel('Overlap')
    plt.savefig(plot_outfile)
    plt.close()

    return 1



def parse_mapping_file(lines, groups_list, delim='\t'):
    labels = lines[0].split(delim)
    indx_sampleid = ''

    group_categ_dict = dict()  # dictionary of dictionary;  the key of dictionary1 = group & value = dictionary2 ;\
                                                           #the key of dictionary2 = category within the group & value = samples belonging to the category within the group

    try:
        indx_sampleid = labels.index("#SampleID")
    except ValueError:
        print "SampleID not in list."

    for g in groups_list:
        indx_categ = ''
        try:
            indx_categ = labels.index(g)
            #print indx_categ
        except ValueError:
            print "%s not in list." %g

        group_categ_dict[g] = list_to_dict(lines, delim, ',', "current", indx_categ, indx_sampleid)

    return group_categ_dict



def main(args):

    parser = argparse.ArgumentParser(description='Walk through the sorted otu lists and plot number of overlapping same taxonomy/OTU')
    parser.add_argument('-i', '--infile1')
    parser.add_argument('-j', '--infile2')
    parser.add_argument('-m', '--mapfile1') # file containing the group where each sample belongs
    parser.add_argument('-n', '--mapfile2') # file containing the group where each sample belongs
    parser.add_argument('-o', '--outfile', default='miseq-JI_count_combo-.png')  # output filename for plot
    parser.add_argument('-s', '--sortedoutfilesuffix', default="-sorted.txt")
    parser.add_argument('-l', '--unlogtransform', type=int) # unlog each value to the base
    parser.add_argument('-t', '--threshold', default=0.001, type=float) # minimum abundance of OTU across the category in a group to be included (default value 0.1% was for css normalized files which are log2 based values) 
    parser.add_argument('-g', '--groups', nargs='*', default=["Treatment"])  # (enter as space delimited on cmd line); creates a list of items used for summarizing
    parser.add_argument('-d', '--delimiter', default='\t')  # delimiter for file
    parser.add_argument('-e', '--extractimpotus', action='store_true', default=False)
    
    args = parser.parse_args()

    if len(sys.argv) <= 4 :
        parser.print_help()
        sys.exit('\natleast five arguments required\n')

    infile1  = read_file(args.infile1)
    infile2  = read_file(args.infile2)

    groups_list = args.groups

    mapfile1  = read_file(args.mapfile1)
    map_file1_dict = parse_mapping_file(mapfile1, groups_list)

    mapfile2  = read_file(args.mapfile2)
    map_file2_dict = parse_mapping_file(mapfile2, groups_list)

    #print map_file1_dict
    #print map_file2_dict

    plot_outfile = args.outfile
    delim = args.delimiter
    suffix = args.sortedoutfilesuffix

    threshold = args.threshold

    log_base = ''
    if args.unlogtransform: # the input data was log transformed?
        log_base = args.unlogtransform

    otu_id_col_file1 = taxonomy_col_file1 = otu_id_col_file2 = taxonomy_col_file2 = '' # ids for the columns

    file1_headers_list = list()
    file2_headers_list = list()

    file1_dict = dict()
    file2_dict = dict()

    file1_otusum_dict = dict()
    file2_otusum_dict = dict()



    for l in infile1:
        l = l.strip('\n')
        if not l:
            continue
        contents = l.split(delim)
        if '#' in l:
            if '#OTU ID' in contents:
                otu_id_col_file1 = contents.index('#OTU ID')
            if 'taxonomy' in contents:
                taxonomy_col_file1 = contents.index('taxonomy')
            file1_headers_list.append(l)
            continue
        
        key = contents[otu_id_col_file1]
        value = ''
        try:
            value= [float(contents[col]) for col in np.arange(otu_id_col_file1+1, taxonomy_col_file1)] # start after otu id and end before taxonomy
        except:
            print l
            sys.exit("A")
        if key not in file1_dict:
            file1_otusum_dict[key] = sum(value)
            value.append(contents[taxonomy_col_file1]) # add taxonomy to the list
            file1_dict[key] = value
        else:
            print l
            sys.exit("B")
        

    for l in infile2:
        l = l.strip('\n')
        if not l:
            continue
        contents = l.split(delim)
        if '#' in l:
            if '#OTU ID' in contents:
                otu_id_col_file2 = contents.index('#OTU ID')
            if 'taxonomy' in contents:
                taxonomy_col_file2 = contents.index('taxonomy')
            file2_headers_list.append(l)
            continue
        
        key = contents[otu_id_col_file2]
        value = ''
        try:
            value= [float(contents[col]) for col in np.arange(otu_id_col_file2+1, taxonomy_col_file2)] # start after otu id and end before taxonomy
        except:
            print l
            sys.exit("C")
        if key not in file2_dict:
            file2_otusum_dict[key] = sum(value)
            value.append(contents[taxonomy_col_file2]) # add taxonomy to the list
            file2_dict[key] = value
        else:
            print l
            sys.exit("D")       

    
    
    #http://stackoverflow.com/questions/33195384/python-sorting-a-dictionary-based-on-values-and-order-based-on-key-if-values-ar

    ListA = list()
    ListB = list()

    sorted_outpfile1 = args.infile1+suffix

    out = open(sorted_outpfile1, 'w')
    for l in file1_headers_list:
        out.write ('%s\n' %(l))

    # for cases where values are same for multiple keys, the keys are asc sorted in lexical order
    for k in sorted(file1_otusum_dict, key=lambda x: (-file1_otusum_dict.get(x), x)):
        ListA.append(k)
        val = [str(j) for j in file1_dict[k]]
        out.write ('%s%s%s\n' %(k, delim, delim.join(val)))
    out.close()

    sorted_outpfile2 = args.infile2+suffix

    out = open(sorted_outpfile2, 'w')
    for l in file2_headers_list:
        out.write ('%s\n' %(l))

    # for cases where values are same for multiple keys, the keys are asc sorted in lexical order
    for k in sorted(file2_otusum_dict, key=lambda x: (-file2_otusum_dict.get(x), x)):
        ListB.append(k)
        val = [str(j) for j in file2_dict[k]]
        out.write ('%s%s%s\n' %(k, delim, delim.join(val)))
    out.close()

    plot_JI_of_files(ListA, ListB, plot_outfile)
    

    top_otus_all_categ_all_groups = []

    output_dumper = list()

    for g , d in map_file1_dict.items():
        #print "File 1" , g
        output_dumper.append("File 1 " + g)
        l = find_top_otus(infile1, d, threshold, delim, log_base)
        output_dumper.extend(l)
        top_otus_all_categ_all_groups.extend(list(set(l)))

    for g , d in map_file2_dict.items():
        #print "File 2" , g
        output_dumper.append("File 2 " + g)
        l = find_top_otus(infile2, d, threshold, delim, log_base)
        output_dumper.extend(l)
        top_otus_all_categ_all_groups.extend(list(set(l)))

    #print list(set(top_otus_all_categ_all_groups))

    # to debug
    #writeLIST_to_file(output_dumper, "tmp-output_dumper.txt")
    #writeLIST_to_file(list(set(top_otus_all_categ_all_groups)), "tmp-top_otus_all_categ_all_groups.txt")

    if args.extractimpotus:
        find_cutoff_number_from_file(set(top_otus_all_categ_all_groups), ListA, ListB)

        extract_important_otus(sorted_outpfile1, set(top_otus_all_categ_all_groups), threshold, '---'.join(groups_list), delim, log_base)
        extract_important_otus(sorted_outpfile2, set(top_otus_all_categ_all_groups), threshold, '---'.join(groups_list), delim, log_base)


def extract_important_otus(ifile, imp_otus, threshold, group_str, delim, log_base):
    infile  = read_file(ifile)

    ofile = ifile+"-groups-"+group_str+"-threshold-"+str(threshold)+"-unloged-lbase-" + str(log_base) +"-imp_otus.txt"
    ofile_orig = ifile+"-groups-"+group_str+"-threshold-"+str(threshold)+ "-orig" +"-imp_otus.txt"

    outlist = list()
    outlist_orig = list()

    otu_id_col_ = taxonomy_col_ = ''

    for l in infile:
        l = l.strip()
        if not l:
            continue

        contents = l.split(delim)  
        if '#' in l:
            outlist.append(l)
            outlist_orig.append(l)
            if '#OTU ID' in contents:
                otu_id_col_ = contents.index('#OTU ID')
            if 'taxonomy' in contents:
                taxonomy_col_ = contents.index('taxonomy')
            continue

        if contents[otu_id_col_] in imp_otus:
            outlist_orig.append(l)
            if log_base:
                # start after otu id and end before taxonomy and convert to log_base^number
                vals = [str(math.pow(log_base,float(contents[col]))) for col in np.arange(otu_id_col_ + 1, taxonomy_col_)] 

                vals.insert(0, contents[otu_id_col_])
                vals.append(contents[taxonomy_col_].replace(":", ";"))

                outlist.append(delim.join(vals))
                
    if log_base: # this file contains values that have been unloged
        writeLIST_to_file(outlist, ofile)
        print "The important OTUs (unlogged values) have been extracted in " + ofile

    # this file contains the original values
    writeLIST_to_file(outlist_orig, ofile_orig)
    print "The important OTUs (original values) have been extracted in " + ofile_orig

    return 1



'''
def find(searchList, elem):
    result_list = []
   
    for e in elem:
        l_list = []
        for i, x in enumerate(searchList):
            if x == e:
                l_list.append(i)
        if len(l_list) == 0: # e was not present in the searchList
            print "%s is not present in the searchList" %e
        result_list.append(l_list)

    return result_list
'''

def find_cutoff_number_from_file(elem, searchListA, searchListB):
    find = lambda searchList, elem: [[i for i, x in enumerate(searchList) if x == e] for e in elem] #http://code.activestate.com/recipes/577943-find-multiple-elements-in-a-list/

    # returns a list of list containing the indexes at which the entries in elem were found in the searchList
    # empty internal list when an entry from elem was not found in searchList
    rankA = find(searchListA, elem)
    rankB = find(searchListB, elem)
    #print elem, rankA, rankB

    #for sublist in l: for item in sublist: yield item
    unpacked_rankA = [item for sublist in rankA for item in sublist]  
    unpacked_rankB = [item for sublist in rankB for item in sublist]

    # since some otuids are unique to each list, the length of elem and unpacked ranks do not match
    #print len(elem), len(searchListA), len(searchListB), len(set(unpacked_rankA)), len(set(unpacked_rankB)) #,min(unpacked_rankA), min(unpacked_rankB)
    print "\nImportant OTUs: %s\nTotal OTUs in file1: %s\nTotal OTUs in file2: %s\nImportant OTUs found in file1: %s\nImportant OTUs found in file2: %s" %(len(elem), len(searchListA), len(searchListB), len(set(unpacked_rankA)), len(set(unpacked_rankB)))

    print "%sTAKE THE TOP %s and %s FROM FILE 1 and 2 TO GET THE TOP OTUs FOR ANY POTENTIAL COMPARISON%sOR" % ('\n', str(max(set(unpacked_rankA))), str(max(set(unpacked_rankB))), '\n')

    print "TAKE THE TOP %s FROM BOTH SORTED FILES TO GET THE TOP OTUs FOR ANY POTENTIAL COMPARISON%s" % (str(max(set(unpacked_rankA + unpacked_rankB))), '\n')


def find_top_otus(lines, categ_samples_dict, threshold, delim, log_base):
   
    # collect samples belonging to the same category within a group
    categ_col_indx_dict = dict() # 0 based indexes, key is category, value is cols of samples in that category

    categ_taxon_values_dict = dict()

    top_otus_all_categ = []


    for k in categ_samples_dict:
        categ_taxon_values_dict[k] = dict()
        taxon_values_dict = dict()
        taxonomy_idx = ''
        sum_otu_abund = 0

        for l in lines:
            l = l.strip()
            contents = l.split(delim)
            if not l or "Constructed from biom file" in l:
                continue
            if '#' in l and '#OTU ID' in l:
                taxonomy_idx = contents.index('#OTU ID')

                # get column numbers of all samples belonging to a category
                v = categ_samples_dict[k]
                categ_col_indx_dict[k] = [contents.index(i) for i in v.split(',') if i in contents]
                continue

            taxon_name = contents[taxonomy_idx]
            v = categ_col_indx_dict[k]

            vals = ''
            if not log_base: # it is empty
                vals = [float(contents[i]) for i in v]
            else:
                vals = [math.pow(log_base,float(contents[i])) for i in v] ## convert to log_base^number

            categ_taxon_values_dict[k][taxon_name] = np.mean(vals)
            taxon_values_dict[taxon_name] = np.mean(vals) 
            sum_otu_abund +=  np.mean(vals)
            #print k, taxon_name, str(np.mean(vals))

        for t in sorted(taxon_values_dict, key=lambda x: (-taxon_values_dict.get(x), x)):
            #print k, t, taxon_values_dict[t]
            top_otus_all_categ.append(t)
            # quit if ratio of otu compared to all otus is less than threshold
            if taxon_values_dict[t]/float(sum_otu_abund) < threshold:
                break
    
    #print categ_taxon_values_dict

    return top_otus_all_categ


###################    


if __name__=='__main__':
    datetime = strftime("%a, %d %b %Y %I:%M:%S %p", localtime())
    cmd = 'echo ' + datetime 
    os.system(cmd)
    main(sys.argv)

