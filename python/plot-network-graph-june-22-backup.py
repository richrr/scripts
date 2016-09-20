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




#python ~/scripts/python/plot-network-graph.py L3_invasive\ 0.310-randomwalk-test.sif L3_invasive\ 0.310\ node_attribute-randomwalk.txt out
#python ~/scripts/python/plot-network-graph.py -i L3_invasive\ 0.310\ edge_attribute-randomwalk.txt -n L3_invasive\ 0.310\ node_attribute-randomwalk.txt
#python ~/scripts/python/plot-network-graph.py -i L3_native\ 0.310\ edge_attribute-randomwalk.txt -n L3_native\ 0.310\ node_attribute-randomwalk.txt

def create_pretty_node_labels(G):
    new_labels = dict()
    for old_label in G.nodes():
        cont = old_label.split("c__")
        new_labels[old_label] = cont[-1] if len(cont) > 1 else cont[0] # UnassignedOtherOther -> cont[0]
        tmp = new_labels[old_label]
        new_labels[old_label] =  tmp if len(tmp)<=5 else tmp[:5]
    return new_labels

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

def identify_node_categ(node_module_dict):
    node_categ_dict = dict()
    for node in node_module_dict:
        start = 'p__'
        end = 'c__'
        result = node
        if start in node and end in node:
            result = re.search('%s(.*)%s' % (start, end), node).group(1)
        node_categ_dict[node] = result
        #print node, result
    return node_categ_dict


def main(args):

    parser = argparse.ArgumentParser(description='Plot networks obtained from MENA')
    parser.add_argument('-i', '--infile') # file containing the edges (interactions between OTUs), can be sif or edge attribute file
    parser.add_argument('-n', '--nodefile') # file containing the node attribute file
    parser.add_argument('-o', '--outfilestr')  # string for output filename
    parser.add_argument('-l', '--logfile', default="ITS-log-file.txt")  # log filename
    parser.add_argument('-q', '--imagequality', default=600, type=int) # image quality dpi
    parser.add_argument('-f', '--imageformat', default='pdf')  # generate images in format. allowed: pdf, png, jpg, tiff
    parser.add_argument('-y', '--outdir', default='./')  # dir name for outputting figs.
    parser.add_argument('-e', '--edgetypeoff', action='store_true', default=False) # do not distinguish positive and negative correlation edges
    parser.add_argument('-d', '--delimiter', default='\t')  # delimiter for file

    args = parser.parse_args()

    if args.infile == None or args.nodefile == None :
        parser.print_help()
        sys.exit('\natleast two arguments (edge file, node attributes file) required\n')

    infile  = args.infile
    node_attrib_file = args.nodefile
    outfilestring = infile.replace(' ','_')
    if args.outfilestr != None:
        outfilestring = args.outfilestr
    outfilestring = args.outdir + outfilestring


    img_qual = args.imagequality
    img_frmt = args.imageformat
    delim = args.delimiter

    #http://www.discoveryplayground.com/computer-programming-for-kids/rgb-colors/
    # purple, deep pink, red, orange, brown, wheat, yellow, forest green, cyan, blue
    attrib_color_map = {'0' : '#a020f0' , '1' : '#ff1493', '2' : '#ff0000', '3' : '#ffa500' , \
                '4' : '#a52a2a' , '5' : '#f5deb3', '6' : '#ffff00' , '7' : '#228b22' , '8' : '#00ffff' , '9' : '#0000ff'}

    #Hot Pink,Violet,Purple,Brown,Salmon,Peru,Orange,Red,Tomato,Blue,Dodger Blue,Deep Sky Blue,Turquoise,Cyan,Light Cyan,Cadet Blue,Mint Cream,Azure,Alice Blue,Lavender,Lavender Blush,Misty Rose,Medium Aquamarine,Aquamarine,Dark Green,Dark Olive Green,Dark Sea Green,Yellow Green,Forest Green,Olive Drab,Dark Khaki,Khaki,Yellow,Light Gray,Green Yellow,Linen,Antique White,Antique White 2,Antique White 3,Antique White 4,Papaya Whip,Blanched Almond,Bisque,Bisque 2,Bisque 3,Bisque 4,Peach Puff,Peach Puff 2,Peach Puff 3,Peach Puff 4,Navajo White,Moccasin,Cornsilk,Cornsilk 2,Cornsilk 3,Cornsilk 4,Ivory,Seashell 2,Seashell 3,Seashell 4,Honeydew,Honeydew 2,Honeydew 3,Honeydew 4,Steel Blue,Light Steel Blue,Light Blue,Powder Blue
    categ_color_map = {'Acidobacteria' : '#ff69b4','Actinobacteria' : '#ee82ee','Aquificae' : '#a020f0',\
    'Armatimonadetes' : '#a52a2a','Bacteroidetes' : '#fa8072','Caldiserica' : '#cd853f','Chlamydiae' : '#ffa500',\
    'Chlorobi' : '#ff0000','Chloroflexi' : '#ff6347','Chrysiogenetes' : '#0000ff','Cyanobacteria' : '#1e90ff',\
    'Deferribacteres' : '#00bfff','Deinococcus-Thermus' : '#40e0d0','Dictyoglomi' : '#00ffff',\
    'Elusimicrobia' : '#e0ffff','Fibrobacteres' : '#5f9ea0','Firmicutes' : '#f5fffa','Fusobacteria' : '#f0ffff',\
    'Gemmatimonadetes' : '#f0f8ff','Lentisphaerae' : '#e6e6fa','Nitrospira' : '#fff0f5','Planctomycetes' : '#ffe4e1',\
    'Proteobacteria' : '#66cdaa','Spirochaetes' : '#7fffd4','Synergistetes' : '#006400','Tenericutes' : '#556b2f',\
    'Thermodesulfobacteria' : '#8fbc8f','Thermomicrobia' : '#9acd32','Thermotogae' : '#228b22','Verrucomicrobia' : '#6b8e23',\
    'Crenarchaeota' : '#bdb76b','Euryarchaeota' : '#f0e68c','Korarchaeota' : '#ffff00','Nanoarchaeota' : '#d3d3d3',\
    'Thaumarchaeota' : '#adff2f','[Parvarchaeota]' : '#faf0e6','[Caldithrix]' : '#faebd7','[Thermi]' : '#eedfcc','AD3' : '#cdc0b0','BHI80-139' : '#8b8378','BRC1' : '#ffefd5','FBP' : '#ffebcd','FCPU426' : '#ffe4c4','GAL15' : '#eed5b7','GN02' : '#cdb79e','GN04' : '#8b7d6b','GOUTA4' : '#ffdab9','Kazan-3B-28' : '#eecbad','MVP-21' : '#cdaf95','MVS-104' : '#8b7765','NC10' : '#ffdead','Nitrospirae' : '#ffe4b5','NKB19' : '#fff8dc','OD1' : '#eee8dc','OP11' : '#cdc8b1','OP3' : '#8b8878','SBR1093' : '#fffff0','SC4' : '#eee5de','SR1' : '#cdc5bf','TM6' : '#8b8682','TM7' : '#f0fff0','WPS-2' : '#e0eee0','WS2' : '#c1cdc1','WS3' : '#838b83','WS4' : '#4682b4','ZB3' : '#b0c4de','Other' : '#add8e6','UnassignedOtherOther' : '#b0e0e6'}


    edge_color_map = {'0':'black', '1.000':'green' , '-1.000':'red'}

    if args.edgetypeoff:
        draw_plots_wout_edge_attributes(infile, node_attrib_file, outfilestring, img_qual, img_frmt, delim, attrib_color_map, categ_color_map, edge_color_map)
    else:
        draw_plots_with_edge_attributes(infile, node_attrib_file, outfilestring, img_qual, img_frmt, delim, attrib_color_map, categ_color_map, edge_color_map)


def get_edge_attributes(edges, edge_color_map):
    edges_attrib_dict = dict()
    if len(edges[0]) < 3:
        sys.exit('\nedge attributes file required, not sif\n')
    else:
        for (u,v,w) in edges:
            edge = (u,v)
            edges_attrib_dict[edge] = edge_color_map[w] # the value is the color based on the edge weight
            #if edge_color_map[w] == 'black':
            #    print edge, '--->', edge_color_map[w]
    return edges_attrib_dict

def draw_plots_with_edge_attributes(infile, node_attrib_file, outfilestring, img_qual, img_frmt, delim, attrib_color_map, categ_color_map, edge_color_map):
    G = nx.Graph()
    # edges
    edges = readColumnsSep(infile, ' ', 0, 2, 4)
    edges_attrib_dict = get_edge_attributes(edges, edge_color_map)
    edges_list = edges_attrib_dict.keys()
    Gcolors = edges_attrib_dict.values()

    #G.add_edges_from(edges_list) # this will only add edges but not the color info. which I might want later
    for (u,v) in edges_list:
        G.add_edge(u , v, color=edges_attrib_dict[(u,v)])

    # node attributes
    node_module_dict = readDict(node_attrib_file, 1, 7, '\t')
    for node in node_module_dict:
        G.add_node(node, moduleNumb = node_module_dict[node])

    #print G.edges(data=True)


    # preprocessed graph, edges only if both nodes are in same module
    edges_between_nodes_from_same_module = preprocess_graph(edges_list, node_module_dict)
    samModG = nx.Graph()

    #samModG.add_edges_from(edges_between_nodes_from_same_module)
    #samModG_edges_list = []
    samModGcolors = []
    for (u,v) in edges_between_nodes_from_same_module:
        # for singleton nodes in module, dummy self edge was added to display on plot
        color_ = ''
        if (u,v) in edges_attrib_dict:
            color_ = edges_attrib_dict[(u,v)]
        elif u==v:
            color_ = 'black'
        samModG.add_edge(u , v, color=color_)
        #samModG_edges_list.append()
        samModGcolors.append(color_)
   
    #print samModGcolors

    # identify the category the OTU belongs to
    node_categ_dict = identify_node_categ(node_module_dict)
    for node in samModG.nodes():
        samModG.add_node(node, category = node_categ_dict[node])

    # reduce length of label for easier visualization
    new_labels = create_pretty_node_labels(G)

    # plot edges as per modules

    #http://stackoverflow.com/questions/24662006/python-networkx-graph-different-colored-nodes-using-two-lists
    nodeColor = [attrib_color_map[G.node[node]['moduleNumb']] for node in G if node != 'Name']

    #http://stackoverflow.com/questions/22992009/legend-in-python-networkx
    # create legend

    f = plt.figure(1)
    ax = f.add_subplot(1,1,1)
    for label in attrib_color_map:
        if label in node_module_dict.values():# only show legend for values that are in my data
                ax.plot([],[],'o',color=attrib_color_map[label],label=label)
    for label in edge_color_map: # 0, 1, -1 correlation values
        if edge_color_map[label] in Gcolors:# colors,  only show legend for values that are in my data
                ax.plot([],[],color=edge_color_map[label],label=label)


    ########### to do ######### calculate relative abundance of the node and make a list and submit as node_size to draw()

    plt.title('OTUs colored as per modules. Intermodule edges allowed.')
    nx.draw(G, edgelist=edges_list, edge_color = Gcolors, node_color = nodeColor, with_labels = False)#, style='dashed')
    #nx.draw_circular(G, node_color = nodeColor, labels = new_labels, with_labels = True)
    #plt.legend()
    plt.legend(bbox_to_anchor=(0.05, 0.93), loc=0, borderaxespad=0.,prop={'size':6}) #, title = "Legend"
    plt.savefig(outfilestring + "-all-edge-node-color-module.png", dpi = img_qual)
    plt.clf()

    # nodecolor as per the phyla
    nodeColor = [categ_color_map[samModG.node[node]['category']] for node in samModG]
    new_labels = create_pretty_node_labels(samModG)

    # create legend
    f = plt.figure(1)
    ax = f.add_subplot(1,1,1)
    for label in categ_color_map:
        if label in node_categ_dict.values():# only show legend for values that are in my data
            ax.plot([],[],'o',color=categ_color_map[label],label=label)
    for label in edge_color_map: # 0, 1, -1 correlation values
        if edge_color_map[label] in samModGcolors:# colors,  only show legend for values that are in my data
                ax.plot([],[],color=edge_color_map[label],label=label)

    plt.title('OTUs colored as per Phylum. Intermodule edges NOT allowed.')
    # other layout algos: dot, neato, fdp, twopi, circo
    algo = 'circo'
    pos = nx.graphviz_layout(samModG, prog=algo)
    nx.draw(samModG, edgelist=edges_between_nodes_from_same_module, edge_color = samModGcolors, pos=pos, node_color = nodeColor, labels = new_labels, with_labels = True)
    #http://stackoverflow.com/questions/7125009/how-to-change-legend-size-with-matplotlib-pyplot 
    #http://matplotlib.org/api/pyplot_api.html#matplotlib.pyplot.legend
    #plt.legend(loc=3,prop={'size':6}) 
    plt.legend(bbox_to_anchor=(0.15, 0.93), loc=0, borderaxespad=0.,prop={'size':6}) #, title = "Legend"
    plt.savefig(outfilestring + "-same-module-edge-node-color-phyla." + img_frmt, dpi = img_qual)
    plt.close()


def draw_plots_wout_edge_attributes(infile, node_attrib_file, outfilestring, img_qual, img_frmt, delim, attrib_color_map, categ_color_map,edge_color_map):
    G = nx.Graph()
    # edges
    edges = readColumnsSep(infile, '\t', 0, 2)
    G.add_edges_from(edges)
    # node attributes
    node_module_dict = readDict(node_attrib_file, 1, 7, '\t')
    for node in node_module_dict:
        G.add_node(node, moduleNumb = node_module_dict[node])

    # preprocessed graph, edges only if both nodes are in same module
    edges_between_nodes_from_same_module = preprocess_graph(edges, node_module_dict)
    samModG = nx.Graph()
    samModG.add_edges_from(edges_between_nodes_from_same_module)
    # identify the category the OTU belongs to
    node_categ_dict = identify_node_categ(node_module_dict)
    for node in samModG.nodes():
        samModG.add_node(node, category = node_categ_dict[node])

    # reduce length of label for easier visualization
    new_labels = create_pretty_node_labels(G)

    # plot edges as per modules

    #http://stackoverflow.com/questions/24662006/python-networkx-graph-different-colored-nodes-using-two-lists
    nodeColor = [attrib_color_map[G.node[node]['moduleNumb']] for node in G if node != 'Name']

    #http://stackoverflow.com/questions/22992009/legend-in-python-networkx
    # create legend

    f = plt.figure(1)
    ax = f.add_subplot(1,1,1)
    for label in attrib_color_map:
        if label in node_module_dict.values():# only show legend for values that are in my data
                ax.plot([0],[0],color=attrib_color_map[label],label=label)

    plt.title('OTUs colored as per modules. Intermodule edges allowed.')
    nx.draw(G, node_color = nodeColor, with_labels = False)
    #nx.draw_circular(G, node_color = nodeColor, labels = new_labels, with_labels = True)
    #plt.legend()
    plt.legend(bbox_to_anchor=(0.05, 0.93), loc=0, borderaxespad=0.,prop={'size':6}) #, title = "Legend"
    plt.savefig(outfilestring + "-all-edge-node-color-module.png", dpi = img_qual)
    plt.clf()

    # nodecolor as per the phyla
    nodeColor = [categ_color_map[samModG.node[node]['category']] for node in samModG]
    new_labels = create_pretty_node_labels(samModG)

    # create legend
    f = plt.figure(1)
    ax = f.add_subplot(1,1,1)
    for label in categ_color_map:
        if label in node_categ_dict.values():# only show legend for values that are in my data
            ax.plot([0],[0],color=categ_color_map[label],label=label)

    plt.title('OTUs colored as per Phylum. Intermodule edges NOT allowed.')
    # other layout algos: dot, neato, fdp, twopi, circo
    algo = 'circo'
    pos = nx.graphviz_layout(samModG, prog=algo)
    nx.draw(samModG, pos=pos, node_color = nodeColor, labels = new_labels, with_labels = True)
    #http://stackoverflow.com/questions/7125009/how-to-change-legend-size-with-matplotlib-pyplot 
    #http://matplotlib.org/api/pyplot_api.html#matplotlib.pyplot.legend
    #plt.legend(loc=3,prop={'size':6}) 
    plt.legend(bbox_to_anchor=(0.15, 0.93), loc=0, borderaxespad=0.,prop={'size':6}) #, title = "Legend"
    plt.savefig(outfilestring + "-same-module-edge-node-color-phyla." + img_frmt, dpi = img_qual)
    plt.close()


if __name__=='__main__':
    datetime = strftime("%a, %d %b %Y %I:%M:%S %p", localtime())
    cmd = 'echo ' + datetime 
    os.system(cmd)
    main(sys.argv)

