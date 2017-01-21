

import os
import sys
from utils import *
import operator
from time import localtime, strftime
import argparse
import re
import numpy as np
import math

# /nfs3/PHARM/Morgun_Lab/richrr/Cervical_Cancer/analysis/single_deg_prob_per_gene
# Usage: python ~/Morgun_Lab/richrr/scripts/python/keep_single_deg_probe_per_gene.py -i ~/Morgun_Lab/richrr/Cervical_Cancer/Data/gene-expression/HumanWG6-IlluminaGeneExpr_QuantileNormalized_Cohort1.csv -g ~/Morgun_Lab/richrr/Cervical_Cancer/analysis/merged/comp/gexpress/p0.2/per_analysis/Analys1-p0.2-consis.csv-comb-pval-output.csv.combpval.0.05.combfdr.0.1.cutoff.csv.consis_genes.csv -m ~/Morgun_Lab/richrr/Cervical_Cancer/Data/gene-expression/common-probes-to-genes-HumanWG6_HumanHT12v4.txt

'''
 -- for each cohort separately(all samples together regardless of the tumor grade):
    - need files:
      -- a data file
      -- list of deg probes
      -- probe to gene symbol mapping
    - among all probes matching the same gene symbol
      -- keep the probe with max expression across all samples
    - output is a list of selected probes (which will be used instead of the consistent gene list)

 -- currently only common ids are kept for down stream analysis. So if the selected (max) probe for
    the same gene is different in the two datasets, the probe will be discarded during merging.
    -- so do the above step for each dataset (cohort) separately, then for non-matching
       selected probes, replace them with the gene symbol
       - need files:
         -- selected probe list from each dataset
         -- data file from each dataset
         -- probe to gene symbol mapping
       - output is list of selected probes (and gene symbols) that are common between two cohorts
         and data files where the different probes for the same gene have been renamed with gene symbol
         in the corresponding data files
'''

def main(args):

    parser = argparse.ArgumentParser(description='among all deg probes matching the same gene symbol keep the deg probe with max expression across all samples')
    parser.add_argument('-i', '--infile')
    parser.add_argument('-g', '--degfile')
    parser.add_argument('-m', '--mapfile')
    
    parser.add_argument('-c', '--colheader', default="PROBE_ID")  # column header when '#' is not in the header line
    parser.add_argument('-o', '--outstr', default="single_deg_probe_")  # output filename
    parser.add_argument('-f', '--fromcolumn', default=0, type=int) # 0-index based column, convert ids FROM this column of mapping file
    parser.add_argument('-t', '--tillcolumn', default=1, type=int) # 0-index based column, convert ids TO this column of mapping file
    parser.add_argument('-r', '--renamecolumn', nargs='*', type =int, default=[0])  # creates a list of items (0-index based column) used as probe ids (this becomes the ids in the FROM column)
    parser.add_argument('-d', '--delimiter', default=',')  # delimiter for file
    parser.add_argument('-s', '--substitute', action='store_true', default=False) # substitute the infile with the new ids from the mapping file 
    
    
    
    args = parser.parse_args()

    if len(sys.argv)==3 :
        parser.print_help()
        sys.exit('\natleast three arguments required\n')

    infile  = args.infile
    degfile = args.degfile
    mapfile = args.mapfile
    
    colheader = args.colheader
    outstr = args.outstr
    fromcol = args.fromcolumn
    tillcol = args.tillcolumn
    renamecol = args.renamecolumn
    delim = args.delimiter

    f2 = read_file(mapfile) # used to make a dictionary of gene symbol (key) to probe id (value)
    d2 = list_to_dict(f2, delim, "__", "current", tillcol, fromcol)

    
    if args.substitute:
        ofile = infile + "_new_id" + ".csv"
        newlist = map_identifiers(infile, mapfile, delim , "__", fromcol, tillcol, False,1, True) # skip missing ids else it throws error in R script if row names have NA 
        #print(newlist)
        writeLIST_to_file(newlist, ofile)
        colheader = newlist[0].split(delim)[0]
        infile = ofile
        print("Changing infile to" + ofile + " and using " + colheader + " as col header")
        

    

    d_list=list()
    comment_list=list()

    for l in open(infile):
        if '#' in l or colheader in l:
            comment_list.append(l.strip())
        else:
            d_list.append(l.strip('\n'))

            


    probe_expression_dict = dict()  # key is probe, value is the sum of expression across all samples
    probe_raw_expression_dict = dict()  # key is probe, value is the expression across all samples
    cannot_convert_to_float = list()
    
    for l in d_list:
        #if "root;" not in l:
        #    continue
        contents = l.strip('\n').split(delim)
        #print contents
        key = delim.join([contents[i] for i in renamecol])
        if key == 'NA':
            key = 'NO_ID_FOUND'
        
        try:
            value= [float(elem) if not 'NA' else 0 for elem in contents[1:] ] # convert all elements except first to float and convert NA to 0
        except:
            #print l
            #sys.exit("Error: Cannot convert to float")
            cannot_convert_to_float.append(l)
            continue
        if key not in probe_raw_expression_dict: # this is not processed so far
            tmp = l.strip('\n').split(delim)[1:] # get everything except first col
            probe_raw_expression_dict[key] = delim.join(tmp)
            probe_expression_dict[key] = sum(value)
        else:
          if args.substitute:
            if sum(value) > probe_expression_dict[key]: # only keep the max valued probe
                tmp = l.strip('\n').split(delim)[1:]
                probe_raw_expression_dict[key] = delim.join(tmp)
                probe_expression_dict[key] = sum(value)
          else:
            print l  # repeated probe
            sys.exit("Error: Repeated probe")


    if args.substitute:
        ofile= infile + outstr + ".csv"
        writeLIST_to_file(comment_list, ofile)  
        writeDICT_to_file(probe_raw_expression_dict, ofile, ',', 'a')
        
        if(len(cannot_convert_to_float) > 0):
            print("\nWarning: Cannot convert these to float so skipped these ids")
            print('\n'.join(cannot_convert_to_float))
        
        sys.exit("\nDone")

   
    degs = [d.strip('\n') for d in read_file(degfile)]
    
    final_list = list()
    # loop over the keys in gene-probe dict and see if elements from deg file
      ## if the gene symbol is '', then use all the deg probes since these novel probes can be researched further
      ## else, add the only probe or probe with max value per gene 
    for k,v in d2.items(): 
        deg_probes = []  # these probes are in the deg file
        for probe in v.split(delim):
            if probe in degs:
                deg_probes.append(probe)

        if len(deg_probes) > 0:
            if not k: # if the gene symbol is empty use all deg probes
                  final_list.extend(deg_probes)
                  print "deg probes with no gene symbol" , len(deg_probes)
            else:
                 selected_deg_probe = deg_probes[0]
                 for p in deg_probes:
                    if probe_expression_dict[p] > probe_expression_dict[ selected_deg_probe ] :
                        selected_deg_probe = p
                 final_list.append( selected_deg_probe )

    outfile = outstr + "consistent.txt"
    writeLIST_to_file(final_list, outfile)
    print "degs after keeping one probe per gene" , len(final_list)
    
    
    



if __name__=='__main__':
    datetime = strftime("%a, %d %b %Y %I:%M:%S %p", localtime())
    cmd = 'echo ' + datetime 
    os.system(cmd)
    main(sys.argv)

