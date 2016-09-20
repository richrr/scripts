import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
from utils import *
import brewer2mpl
import re
import os
import sys
import operator
from time import localtime, strftime
import argparse
import os.path
from subprocess import Popen, PIPE
from collections import defaultdict
import math
from matplotlib_venn import venn2, venn2_circles


'''
# compare the edges between networks
    * present/absence of edges
    * different type of correlation/interaction (1 vs -1)
    * np vs pp vs ...

# compare the module assignment of the node
    * are the modules numbered from max. nodes to min. nodes or vice-versa
    * for each node: scatter plot of module to which it belongs
'''

# python ~/scripts/python/compare-network-graphs.py -i L3_invasive\ 0.310\ edge_attribute-leadeigenvect.txt L3_native\ 0.310\ edge_attribute-leadeigenvect.txt -n L3_invasive\ 0.310\ node_attribute-leadeigenvect.txt L3_native\ 0.310\ node_attribute-leadeigenvect.txt -l L3_inv_nat_0.310_leadeigenvect-log.txt

def preprocess_graph(edges, node_module_dict):
    edges_between_nodes_from_same_module = list()
    for (u,v) in edges:
        #print u, "---->", v
        if u in node_module_dict and v in node_module_dict:
            if node_module_dict[u] == node_module_dict[v]:
                edges_between_nodes_from_same_module.append((u,v))
        else:
            print "One or both nodes of the edge are missing from the attribute file"
    # for nodes that do not have edges when inter module edges are removed
    for (u,v) in edges:
        u_present = [tup for tup in edges_between_nodes_from_same_module if u in tup]
        if len(u_present) == 0:
            edges_between_nodes_from_same_module.append((u,u))
        v_present = [tup for tup in edges_between_nodes_from_same_module if v in tup]
        if len(v_present) == 0:
            edges_between_nodes_from_same_module.append((v,v))
    return edges_between_nodes_from_same_module

def main(args):

    parser = argparse.ArgumentParser(description='Plot networks obtained from MENA')
    parser.add_argument('-i', '--infiles', nargs=2) # files containing the edges (interactions between OTUs), can be sif or edge attribute file
    parser.add_argument('-n', '--nodefiles', nargs=2) # files containing the node attribute
    parser.add_argument('-a', '--nodeabund') # file containing the node abundance
    parser.add_argument('-o', '--outfilestr')  # string for output filename
    parser.add_argument('-l', '--logfile', default="log-file.txt")  # log filename
    parser.add_argument('-q', '--imagequality', default=600, type=int) # image quality dpi
    parser.add_argument('-f', '--imageformat', default='pdf')  # generate images in format. allowed: pdf, png, jpg, tiff
    parser.add_argument('-y', '--outdir', default='./')  # dir name for outputting figs.
    parser.add_argument('-e', '--edgetypeoff', action='store_true', default=False) # do not distinguish positive and negative correlation edges
    parser.add_argument('-d', '--delimiter', default='\t')  # delimiter for file
    parser.add_argument('-t', '--twomethods', action='store_true', default=False) # common edges of same group from two methods
    parser.add_argument('-c', '--comparetwomethodsnets', action='store_true', default=False) # compare same module common edges (nets) of differents group from two methods. P.S. you have to run the script twice with -t and then finally once with -c

    args = parser.parse_args()
    if args.infiles == None:
        parser.print_help()
        sys.exit('\natleast two arguments (edge file, node attributes file or -c) required\n')
    
    if (not args.comparetwomethodsnets) and args.nodefiles == None :
        parser.print_help()
        sys.exit('\nEither node attributes files Or -c argument required\n')

    infiles  = args.infiles

    log_file = args.logfile
    log_file_list = [str(args)]

    img_qual = args.imagequality
    img_frmt = args.imageformat
    delim = args.delimiter

    node_attrib_files = ''
    if args.nodefiles != None:
        node_attrib_files = args.nodefiles

    if args.edgetypeoff:
        sys.exit("No method found\n")
    else:
        comp_networks_with_edge_attributes(infiles, node_attrib_files, img_qual, img_frmt, delim, log_file_list, log_file, args.twomethods, args.comparetwomethodsnets)

def get_edge_correlations(edges):
    edges_correl_dict = dict()
    if len(edges[0]) < 3:
        sys.exit('\nedge attributes file required, not sif\n')
    else:
        for (u,v,w) in edges:
            edge = (u,v)
            edges_correl_dict[edge] = w # the value is the edge weight/correlation
    return edges_correl_dict

def dict_diff(d_first, d_second):
    sd1 = set(d_first)
    sd2 = set(d_second)
    #Edges missing in the second dict
    edgeNOTin2 = sd1.difference(sd2)
    #Edges missing in the first dict
    edgeNOTin1 = sd2.difference(sd1)
    #Check for differences
    common_edges = sd1.intersection(sd2)  # common edges
    common_edges_w_diff_correl = list() # common edges with different correlations
    for edge in common_edges:
        if d_first[edge] != d_second[edge]:
            common_edges_w_diff_correl.append(edge)
    # total edges
    all_edges = sd1.union(sd2)
    return edgeNOTin1, edgeNOTin2, common_edges, all_edges, common_edges_w_diff_correl

def print_network_comparisons(edges_correl_dict1, edges_correl_dict2, edgeNOTin1, edgeNOTin2, common_edges, all_edges, common_edges_w_diff_correl):
    edges_nets = "\nEdges in networks 1, 2: %s, %s" %(len(edges_correl_dict1), len(edges_correl_dict2))
    edges_intersect = "Intersection-> Edges present in both networks: %s" %  (len(common_edges))
    edges_union = "Union-> Edges present in networks 1 OR 2: %s" %(len(all_edges))

    prop_edges_both_nets = "Proportion of edges present in both networks: %s" %  (len(common_edges)/float(len(all_edges)))
    names_edges_both_nets = "Edges present in both networks:\n\t" + '\n\t'.join(["%s,%s" % tup for tup in common_edges])

    edges_both_nets_diff_correl = "Edges present in both networks but have different correlations: %s"   % (len(common_edges_w_diff_correl))
    prop_edges_both_nets_diff_correl =  "Proportion of edges present in both networks but have different correlations: %s"  \
                            % (len(common_edges_w_diff_correl)/float(len(all_edges)))
    names_edges_both_nets_diff_correl = "Edges present in both networks but have different correlations:\n\t" + '\n\t'.join(["%s,%s" % tup for tup in common_edges_w_diff_correl])

    edges_in_1_not_2 = "Edges present in first network but not other: " + str(len(edgeNOTin2)) + "\n" + '\n\t'.join(["%s,%s" % tup for tup in edgeNOTin2])
    edges_in_2_not_1 = "Edges present in second network but not other: " + str(len(edgeNOTin1)) + "\n" + '\n\t'.join(["%s,%s" % tup for tup in edgeNOTin1])

    res = [edges_nets, edges_intersect, edges_union, prop_edges_both_nets, names_edges_both_nets, edges_both_nets_diff_correl, \
			prop_edges_both_nets_diff_correl, names_edges_both_nets_diff_correl, edges_in_1_not_2, edges_in_2_not_1]
    return res

def plot_venn2C(infile1, infile2, lst1, lst2, info='', outfile="venn", img_frmt=".png", img_qual=600):
    f1 = infile1
    f2 = infile2
    old_new_patt = {' edge_attribute':'', '.txt':'', ' ':'\n' }
    for old, new in old_new_patt.items():
        f1 = f1.replace(old, new)
        f2 = f2.replace(old, new)
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.text(0.8, 0.99, info, style='italic', fontsize=10, bbox={'facecolor':'blue', 'alpha':0.5, 'pad':10}, horizontalalignment='right', verticalalignment='top', transform=ax.transAxes)
    #v = venn2([set(lst1), set(lst2)], set_labels = ('f1', 'f2'))
    # the following allows to adjust font size and text of set labels
    v = venn2([set(lst1), set(lst2)])#, set_labels = ('f1', 'f2'))
    A = v.get_label_by_id('A')
    A.set_text(f1)
    A.set_fontsize(10)
    #A.set_x(A.get_position()[0] + 0.1)
    #A.set_y(A.get_position()[1] + 1.0)

    B = v.get_label_by_id('B')
    B.set_text(f2)
    B.set_fontsize(10)
    #B.set_x(B.get_position()[0] - 0.1)
    #B.set_y(B.get_position()[1] + 1.0)
    
    #labelA = v.get_label_by_id('10')          # Those are subset labels (i.e. numbers)
    #print labelA.get_text()
    #labelA.set_text()
    #labelA.set_fontsize(22) 
    #labelA.set_family('serif')
    #labelA.set_x(label.get_position()[0] + 0.1)

    outfilename = f1+"_"+f2
    outfilename = outfilename.replace('\n','_') + '_' + outfile
    outfilename = outfilename.replace('../','')
    outfilename = outfilename.replace('/','_')
    outfilename += img_frmt
    ax.axis()
    plt.savefig(outfilename , dpi = img_qual)
    plt.clf()
    plt.close()
    return outfilename

def comp_networks_with_edge_attributes(infiles, node_attrib_files, img_qual, img_frmt, delim, log_file_list, log_file, twomethods, comparetwomethodsnets):
    infile1 = infile2 = ''
    if len(infiles) == 2:
        infile1 = infiles[0]
        infile2 = infiles[1]

    # edges
    edges1 = readColumnsSep(infile1, ' ', 0, 2, 4)
    edges2 = readColumnsSep(infile2, ' ', 0, 2, 4)
    edges_correl_dict1 = get_edge_correlations(edges1)
    edges_correl_dict2 = get_edge_correlations(edges2)
    # compare two edge lists for presene absence of list:
    edgeNOTin1, edgeNOTin2, common_edges, all_edges, common_edges_w_diff_correl = dict_diff(edges_correl_dict1, edges_correl_dict2)

    res = print_network_comparisons(edges_correl_dict1, edges_correl_dict2, edgeNOTin1, edgeNOTin2, common_edges, all_edges, common_edges_w_diff_correl)
    log_file_list.extend(res)
    
    edges_list1 = edges_correl_dict1.keys()
    edges_list2 = edges_correl_dict2.keys()
    info = 'Edges in both "entire" networks but diff correlations: %s'   % (len(common_edges_w_diff_correl))
    fign = plot_venn2C(infile1, infile2, edges_list1, edges_list2, info)
    add_to_list = 'Venn Fig. saved as %s'%fign
    log_file_list.append(add_to_list)

    fign = fign.replace('.png' , '-') # remove png to reduce file name length

    # write to file: (i) common edges (ii) common edges with diff. correlations
    writeLIST_to_file(["%s,%s" % tup for tup in common_edges], fign+"-commEdges.txt")
    writeLIST_to_file(["%s,%s" % tup for tup in common_edges_w_diff_correl], fign+"-commEdges_w_diff_correl.txt")

    if comparetwomethodsnets:  # run code below when -n was not but -c was provided, directly go to plot networks
        run_cmds = '\nRun two commands: \n'
        for (e, n) in [(infile1, "node_attrib_file1") , (infile2, "node_attrib_file2")]:
            run_cmds += 'python ~/scripts/python/plot-network-graph.py -i %s -n %s -a ?_node_abundance.txt -c ?-commEdges.txt -b ?-commEdges_w_diff_correl.txt\n' %(e, n)
        print run_cmds
        log_file_list.append(run_cmds)
    else: # run code below when -c was not but -n was provided
        # node attributes
        node_attrib_file1 = node_attrib_file2 = ''
        if len(node_attrib_files) == 2:
            node_attrib_file1 = node_attrib_files[0]
            node_attrib_file2 = node_attrib_files[1]
        node_module_dict1 = readDict(node_attrib_file1, 1, 7, '\t')
        node_module_dict2 = readDict(node_attrib_file2, 1, 7, '\t')
    
        # preprocessed graph, edges only if both nodes are in same module 
        edges_between_nodes_from_same_module1 = preprocess_graph(edges_list1, node_module_dict1)
        edges_between_nodes_from_same_module2 = preprocess_graph(edges_list2, node_module_dict2)
        
        # dict of edges between nodes of same module
        sm_edges_dict1 = dict()
        sm_edges_dict2 = dict()
        for e in edges_between_nodes_from_same_module1:
            sm_edges_dict1[e] = edges_correl_dict1[e] if e in edges_correl_dict1 else '0' # to avoid key error for dummy self edges
        
        for e in edges_between_nodes_from_same_module2:
            sm_edges_dict2[e] = edges_correl_dict2[e] if e in edges_correl_dict2 else '0' # to avoid key error for dummy self edges
        
        # calculate above differences for preprocessed networks (where intermodule edges are removed)
        # compare two edge lists for presene absence of list:
        sm_edgeNOTin1, sm_edgeNOTin2, sm_common_edges, sm_all_edges, sm_common_edges_w_diff_correl = dict_diff(sm_edges_dict1, sm_edges_dict2)
    
        res = print_network_comparisons(sm_edges_dict1, sm_edges_dict2, sm_edgeNOTin1, sm_edgeNOTin2, sm_common_edges, sm_all_edges, sm_common_edges_w_diff_correl)
        log_file_list.extend(res)
    
        sm_edges_list1 = sm_edges_dict1.keys()
        sm_edges_list2 = sm_edges_dict2.keys()
        info = 'Edges in both "same module" networks but diff correlations: %s'   % (len(sm_common_edges_w_diff_correl))
        fign = plot_venn2C(infile1, infile2, sm_edges_list1, sm_edges_list2, info, 'venn-sm')
        add_to_list = 'Venn-SM Fig. saved as %s'%fign
        log_file_list.append(add_to_list)

        fign = fign.replace('.png' , '-') # remove png to reduce file name length

        # write to file: (i) common edges (ii) common edges with diff. correlations
        writeLIST_to_file(["%s,%s" % tup for tup in sm_common_edges], fign+"-SM_commEdges.txt")
        writeLIST_to_file(["%s,%s" % tup for tup in sm_common_edges_w_diff_correl], fign+"-SM_commEdges_w_diff_correl.txt")
    
        # Write the two-methods common edges (w and w/o) same correlations, with edge info. 
        # (1 for +ve, -1 for -ve, 99 for diff. correlations)
        if twomethods:
            two_method_edges = list()
            for tup in sm_common_edges:
                edge_info = ''
                if tup not in sm_common_edges_w_diff_correl:
                    edge_info = "%s (na) %s = %s" %(tup[0], tup[1], sm_edges_dict1[tup])
                else:
                    edge_info = "%s (na) %s = 99" %tup
                two_method_edges.append(edge_info)
            writeLIST_to_file(two_method_edges, fign+"-SM_commEdges-2methEdges.txt")
        else: # directly go to plot networks
            run_cmds = '\nRun two commands: \n'
            for (e, n) in [(infile1, node_attrib_file1) , (infile2, node_attrib_file2)]:
                run_cmds += 'python ~/scripts/python/plot-network-graph.py -i %s -n %s -a ?_node_abundance.txt -c ?-SM_commEdges.txt -b ?-SM_commEdges_w_diff_correl.txt\n' %(e, n)
            print run_cmds
            log_file_list.append(run_cmds)

    
    # write to log file
    writeLIST_to_file(log_file_list, log_file)

if __name__=='__main__':
    datetime = strftime("%a, %d %b %Y %I:%M:%S %p", localtime())
    cmd = 'echo ' + datetime 
    os.system(cmd)
    main(sys.argv)

