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
import matplotlib.font_manager as font_manager
import matplotlib


#python ~/scripts/python/plot-network-graph.py -i L3_invasive\ 0.310\ edge_attribute-leadeigenvect.txt -n L3_invasive\ 0.310\ node_attribute-leadeigenvect.txt -a L3_invasive.txt -c L3_invasive_0.310-leadeigenvect_L3_native_0.310-leadeigenvect_venn-sm.png-SM_common_edges.txt -b L3_invasive_0.310-leadeigenvect_L3_native_0.310-leadeigenvect_venn-sm.png-SM_common_edges_w_diff_correl.txt


###### to do: input 4 files: 1 common edge file, 1 common edge file with diff correl, 1 common edge file SM, 1 common edge file with diff correl SM
######## or just input the 1 common edge file SM, 1 common edge file with diff correl SM (since we wouldn't show the edges in the entire network)

def create_pretty_node_labels(G):
    new_labels = dict()
    for old_label in G.nodes():
        pattern_to_find = "c__"
        if 'o__' in old_label:
            pattern_to_find = 'o__'
        cont = old_label.split(pattern_to_find)
        new_labels[old_label] = cont[-1] if len(cont) > 1 else cont[0] # UnassignedOtherOther -> cont[0]
        tmp = new_labels[old_label]
        #new_labels[old_label] =  tmp if len(tmp)<=5 else tmp[:5]
        new_labels[old_label] =  tmp if len(tmp)<=5 else insert_newlines(tmp)
    return new_labels

def insert_newlines(string, every=6):
    return '\n'.join(string[i:i+every] for i in xrange(0, len(string), every))

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
        elif start in node and end not in node:
            result = node[node.index(start)+3:]
        elif start not in node and end not in node:
            if 'k__Fungi' in node:
                result = node[node.index('k__Fungi')+8:] # hack for k__FungiOtherOtherOther
        node_categ_dict[node] = result
        #print node, result
    return node_categ_dict

def create_color_map(node_categ_dict_values):
    categ_color_map = dict()
    available_colors_key = ['#ff69b4','#ee82ee','#a020f0','#a52a2a','#fa8072','#cd853f','#ffa500','#ff0000',\
'#ff6347','#0000ff','#1e90ff','#00bfff','#40e0d0','#00ffff','#e0ffff','#5f9ea0','#f5fffa',\
'#f0ffff','#f0f8ff','#e6e6fa','#fff0f5','#ffe4e1','#66cdaa','#7fffd4','#6400','#556b2f','#8fbc8f',\
'#9acd32','#228b22','#6b8e23','#bdb76b','#f0e68c','#ffff00','#d3d3d3','#adff2f','#faf0e6','#faebd7',\
'#eedfcc','#cdc0b0','#8b8378','#ffefd5','#ffebcd','#ffe4c4','#eed5b7','#cdb79e','#8b7d6b','#ffdab9',\
'#eecbad','#cdaf95','#8b7765','#ffdead','#ffe4b5','#fff8dc','#eee8dc','#cdc8b1','#8b8878','#fffff0',\
'#eee5de','#cdc5bf','#8b8682','#f0fff0','#e0eee0','#c1cdc1','#838b83','#4682b4','#b0c4de','#add8e6','#b0e0e6']
    idx = 0
    for phyla in set(node_categ_dict_values):
        if idx <= len(available_colors_key):
            pass
        else:
            idx = 0
        categ_color_map[phyla] = available_colors_key[idx]
        idx += 1
    return categ_color_map

def main(args):

    parser = argparse.ArgumentParser(description='Plot networks obtained from MENA')
    parser.add_argument('-i', '--infile') # file containing the edges (interactions between OTUs), can be sif or edge attribute file
    parser.add_argument('-n', '--nodefile') # file containing the node attribute
    parser.add_argument('-c', '--commonedgefile') # file containing the edges common with the "other" network (here, used only for the Same Module network)
    parser.add_argument('-b', '--commonedgediffcorrelfile') # file containing the edges common with the "other" network but different correlations/edge attributes.  (here, used only for the Same Module network)
    parser.add_argument('-a', '--nodeabund') # file containing the node abundance
    parser.add_argument('-o', '--outfilestr')  # string for output filename
    parser.add_argument('-l', '--logfile', default="ITS-log-file.txt")  # log filename
    parser.add_argument('-q', '--imagequality', default=600, type=int) # image quality dpi
    parser.add_argument('-f', '--imageformat', default='pdf')  # generate images in format. allowed: pdf, png, jpg, tiff
    parser.add_argument('-y', '--outdir', default='./')  # dir name for outputting figs.
    parser.add_argument('-e', '--edgetypeoff', action='store_true', default=False) # do not distinguish positive and negative correlation edges
    parser.add_argument('-u', '--fungal', action='store_true', default=False) # use the fungal categ_color_map instead of (default) bacteria
    parser.add_argument('-m', '--dynamic', action='store_true', default=False) # use the dynamic categ_color_map instead of (default) bacteria. useful for bacteria-fungi networks.
    parser.add_argument('-k', '--bactfungal', action='store_true', default=False) # use the bactfungal categ_color_map instead of (default) bacteria. useful for bacteria-fungi networks. This is hack so the native and invasive use the same colorscheme and easy to compare.
    parser.add_argument('-d', '--delimiter', default='\t')  # delimiter for file
    parser.add_argument('-z', '--comparetwomethodsnets', action='store_true', default=False) # compare same module common edges (nets) of differents group from two methods. P.S. you have to run the script twice with -t and then finally once with -c

    args = parser.parse_args()

    if args.infile == None:
        parser.print_help()
        sys.exit('\natleast two arguments (edge file, node attributes file or -z) required\n')

    infile  = args.infile

    node_attrib_file = ''
    if args.nodefile != None:
        node_attrib_file = args.nodefile

    outfilestring = infile.replace(' ','_')
    if args.outfilestr != None:
        outfilestring = args.outfilestr
    outfilestring = args.outdir + outfilestring + '_pubs_'
    if args.commonedgefile != None and args.commonedgediffcorrelfile != None:
        outfilestring += 'withCommonEdgeInfo' 
    node_abund_file = args.nodeabund

    #common_edge_list = []
    #common_edge_diff_correl_list = []
    common_edge_list = [i.replace('\n', '') for i in read_file(args.commonedgefile)] if args.commonedgefile != None else []
    common_edge_diff_correl_list = [i.replace('\n', '') for i in read_file(args.commonedgediffcorrelfile)] if args.commonedgediffcorrelfile != None else []

    img_qual = args.imagequality
    img_frmt = args.imageformat
    delim = args.delimiter

    #http://www.discoveryplayground.com/computer-programming-for-kids/rgb-colors/
    # purple, deep pink, red, orange, brown, wheat, yellow, forest green, cyan, blue
    attrib_color_map = {'0' : '#a020f0' , '1' : '#ff1493', '2' : '#ff0000', '3' : '#ffa500' , \
                '4' : '#a52a2a' , '5' : '#f5deb3', '6' : '#ffff00' , '7' : '#228b22' , '8' : '#00ffff' , '9' : '#0000ff'}

    categ_color_map = {}
    if args.dynamic:
        pass
    elif args.bactfungal:
        categ_color_map = {'OtherOtherOther' : '#ff1493','Ascomycota' : '#f5deb3','AscomycotaOtherOther' : '#ffa07a',\
    'Basidiomycota' : '#ff8c00','BasidiomycotaOtherOther' : '#ff4500', 'Chytridiomycota' : '#b22222','unidentified' : '#d2691e',\
    'Zygomycota' : '#ffd700','Acidobacteria' : '#ff69b4','Actinobacteria' : '#ee82ee','Aquificae' : '#a020f0',\
    'Armatimonadetes' : '#a52a2a','Bacteroidetes' : '#fa8072','Caldiserica' : '#cd853f','Chlamydiae' : '#ffa500',\
    'Chlorobi' : '#ff0000','Chloroflexi' : '#ff6347','Chrysiogenetes' : '#0000ff','Cyanobacteria' : '#1e90ff',\
    'Deferribacteres' : '#00bfff','Deinococcus-Thermus' : '#40e0d0','Dictyoglomi' : '#00ffff',\
    'Elusimicrobia' : '#e0ffff','Fibrobacteres' : '#5f9ea0','Firmicutes' : '#f5fffa','Fusobacteria' : '#f0ffff',\
    'Gemmatimonadetes' : '#f0f8ff','Lentisphaerae' : '#e6e6fa','Nitrospira' : '#fff0f5','Planctomycetes' : '#ffe4e1',\
    'Proteobacteria' : '#66cdaa','Spirochaetes' : '#7fffd4','Synergistetes' : '#006400','Tenericutes' : '#556b2f',\
    'Thermodesulfobacteria' : '#8fbc8f','Thermomicrobia' : '#9acd32','Thermotogae' : '#228b22','Verrucomicrobia' : '#6b8e23',\
    'Crenarchaeota' : '#bdb76b','Euryarchaeota' : '#f0e68c','Korarchaeota' : '#ffff00','Nanoarchaeota' : '#d3d3d3',\
    'Thaumarchaeota' : '#adff2f','[Parvarchaeota]' : '#faf0e6','[Caldithrix]' : '#faebd7','[Thermi]' : '#eedfcc','AD3' : '#cdc0b0','BHI80-139' : '#8b8378','BRC1' : '#ffefd5','FBP' : '#ffebcd','FCPU426' : '#ffe4c4','GAL15' : '#eed5b7','GN02' : '#cdb79e','GN04' : '#8b7d6b','GOUTA4' : '#ffdab9','Kazan-3B-28' : '#eecbad','MVP-21' : '#cdaf95','MVS-104' : '#8b7765','NC10' : '#ffdead','Nitrospirae' : '#ffe4b5','NKB19' : '#fff8dc','OD1' : '#eee8dc','OP11' : '#cdc8b1','OP3' : '#8b8878','SBR1093' : '#fffff0','SC4' : '#eee5de','SR1' : '#cdc5bf','TM6' : '#8b8682','TM7' : '#f0fff0','WPS-2' : '#e0eee0','WS2' : '#c1cdc1','WS3' : '#838b83','WS4' : '#4682b4','ZB3' : '#b0c4de','Other' : '#add8e6','UnassignedOtherOther' : '#b0e0e6'}
    elif args.fungal:
        categ_color_map = {'OtherOtherOther' : '#ff69b4','Ascomycota' : '#ee82ee','AscomycotaOtherOther' : '#a020f0',\
    'Basidiomycota' : '#a52a2a','BasidiomycotaOtherOther' : '#fa8072', 'Chytridiomycota' : '#0000ff','unidentified' : '#1e90ff',\
    'Zygomycota' : '#00bfff'}
    else:
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
    #Hot Pink,Violet,Purple,Brown,Salmon,Peru,Orange,Red,Tomato,Blue,Dodger Blue,Deep Sky Blue,Turquoise,Cyan,Light Cyan,Cadet Blue,Mint Cream,Azure,Alice Blue,Lavender,Lavender Blush,Misty Rose,Medium Aquamarine,Aquamarine,Dark Green,Dark Olive Green,Dark Sea Green,Yellow Green,Forest Green,Olive Drab,Dark Khaki,Khaki,Yellow,Light Gray,Green Yellow,Linen,Antique White,Antique White 2,Antique White 3,Antique White 4,Papaya Whip,Blanched Almond,Bisque,Bisque 2,Bisque 3,Bisque 4,Peach Puff,Peach Puff 2,Peach Puff 3,Peach Puff 4,Navajo White,Moccasin,Cornsilk,Cornsilk 2,Cornsilk 3,Cornsilk 4,Ivory,Seashell 2,Seashell 3,Seashell 4,Honeydew,Honeydew 2,Honeydew 3,Honeydew 4,Steel Blue,Light Steel Blue,Light Blue,Powder Blue

   

    edge_color_map = {'0':'cyan', '1.000':'green' , '-1.000':'red' , '99':'commEdgeDiffCorrel2methods', 'commEdge':'black' , 'commEdgeDiffCorrel': 'blue'}

    if args.edgetypeoff:
        draw_plots_wout_edge_attributes(infile, node_attrib_file, outfilestring, img_qual, img_frmt, delim, attrib_color_map, categ_color_map, edge_color_map)
    if args.comparetwomethodsnets:
        draw_plots_with_edge_attributes_no_node_module_file(infile, node_attrib_file, outfilestring, img_qual, img_frmt, delim, attrib_color_map, categ_color_map, edge_color_map, node_abund_file, common_edge_list, common_edge_diff_correl_list, args.dynamic)
    else:
        draw_plots_with_edge_attributes(infile, node_attrib_file, outfilestring, img_qual, img_frmt, delim, attrib_color_map, categ_color_map, edge_color_map, node_abund_file, common_edge_list, common_edge_diff_correl_list)


def get_edge_attributes(edges, edge_color_map, common_edge_list, common_edge_diff_correl_list):
    edges_attrib_dict = dict()
    if len(edges[0]) < 3:
        sys.exit('\nedge attributes file required, not sif\n')
    else:
        for (u,v,w) in edges:
            edge = (u,v)
            if ','.join(edge) in common_edge_diff_correl_list:
                #something
                edges_attrib_dict[edge] = edge_color_map['commEdgeDiffCorrel']
                continue
            if ','.join(edge) in common_edge_list:
                #something
                edges_attrib_dict[edge] = edge_color_map['commEdge']
                continue
            edges_attrib_dict[edge] = edge_color_map[w] # the value is the color based on the edge weight
    return edges_attrib_dict

def get_node_abundance(node_abund_file, delim):
    lines = read_file(node_abund_file)
    node_avgAbund_dict = dict()
    for l in lines:
        if '#' in l or 'Taxon' in l:
            continue
        l = l.strip().split(delim)
        taxon = l[0]
        abund = [float(i) for i in l[1:]]
        avg_abund = sum(abund)/float(len(abund))
        node_avgAbund_dict[taxon] = avg_abund
    return node_avgAbund_dict

def node_weights_to_sizes(x):
    sorted_x = sorted(x.items(), key=operator.itemgetter(1))
    node_sizes_dict = dict()
    for idx, val in enumerate(sorted_x):
        node = val[0]
        for patt in [';', ' ', '[', ']']:
            node = node.replace(patt, '')
        node_sizes_dict[node] = 750 + (idx * 2) # temp hack to increase size of nodes for publication
        #node_sizes_dict[node] = 150 + (idx * 2)
    #print node_sizes_dict
    return node_sizes_dict

def draw_plots_with_edge_attributes_no_node_module_file(infile, node_attrib_file, outfilestring, img_qual, img_frmt, delim, attrib_color_map, categ_color_map, edge_color_map, node_abund_file, common_edge_list, common_edge_diff_correl_list, dynamic):
    # edges
    edges = readColumnsSep(infile, ' ', 0, 2, 4)
    edges_attrib_dict = get_edge_attributes(edges, edge_color_map, common_edge_list, common_edge_diff_correl_list)
    edges_list = edges_attrib_dict.keys()

    # calculate avg. relative abundance of the node
    node_abundance = get_node_abundance(node_abund_file, delim)
    #and make a list and submit as node_size to draw()
    node_sizes_dict = node_weights_to_sizes(node_abundance)

    samModG = nx.Graph()
    
    samModGcolors = []
    for (u,v) in edges_list:
        # for singleton nodes in module, dummy self edge was added to display on plot
        color_ = ''
        if (u,v) in edges_attrib_dict:
            color_ = edges_attrib_dict[(u,v)]
        elif u==v:
            color_ = 'cyan'
        samModG.add_edge(u , v, color=color_)
        samModGcolors.append(color_)
   

    all_nodes_in_edge_list = [','.join(e) for e in edges_list]
    # identify the category the OTU belongs to
    node_categ_dict = identify_node_categ(condense_list(all_nodes_in_edge_list , ','))
    # use to dynamically create color map
    if dynamic:
        categ_color_map = create_color_map(node_categ_dict.values()) 
    for node in samModG.nodes():
        samModG.add_node(node, category = node_categ_dict[node])

    samModGnode_sizes_list = [node_sizes_dict[i] for i in samModG.nodes()]

    # reduce length of label for easier visualization
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

    plt.title('OTUs colored as per Phylum.')
    # other layout algos: dot, neato, fdp, twopi, circo
    algo = 'circo'
    pos = nx.graphviz_layout(samModG, prog=algo)


    #https://wiki.ubuntu.com/Fonts
    #williamslab@HORT-MW-Vos:/usr/share/matplotlib/mpl-data/fonts$ sudo ln -s /usr/share/fonts/truetype/msttcorefonts/Arial.ttf ./ttf/arial.ttf
    #williamslab@HORT-MW-Vos:/usr/share/matplotlib/mpl-data/fonts$ sudo ln -s /usr/share/fonts/truetype/msttcorefonts/Times_New_Roman.ttf ./ttf/times.ttf 
    
    '''
    fontpath = '/usr/local/share/fonts/Arial.ttf'
    prop = font_manager.FontProperties(fname=fontpath)
    matplotlib.rcParams['font.family'] = prop.get_name()
    '''
    #/usr/share/fonts/truetype/msttcorefonts/
    text_font = 'Arial' # 'Times' 'Helvetica'
    '''
    #http://stackoverflow.com/questions/18821795/how-can-i-get-list-of-font-familyor-name-of-font-in-matplotlib
    '''


    nx.draw(samModG, edgelist=edges_list, edge_color = samModGcolors, pos=pos, node_color = nodeColor, labels = new_labels, with_labels = True, node_size=samModGnode_sizes_list, font_size=8, font_family=text_font)
    #http://stackoverflow.com/questions/7125009/how-to-change-legend-size-with-matplotlib-pyplot 
    #http://matplotlib.org/api/pyplot_api.html#matplotlib.pyplot.legend
    #plt.legend(loc=3,prop={'size':6}) 
    plt.legend(bbox_to_anchor=(0.15, 0.93), loc=0, borderaxespad=0.,prop={'size':6}) #, title = "Legend"
    plt.savefig(outfilestring + "-edge-node-color-phyla." + img_frmt, dpi = img_qual)
    plt.close()


def draw_plots_with_edge_attributes(infile, node_attrib_file, outfilestring, img_qual, img_frmt, delim, attrib_color_map, categ_color_map, edge_color_map, node_abund_file, common_edge_list, common_edge_diff_correl_list):
    G = nx.Graph()
    # edges
    edges = readColumnsSep(infile, ' ', 0, 2, 4)
    edges_attrib_dict = get_edge_attributes(edges, edge_color_map, common_edge_list, common_edge_diff_correl_list)
    edges_list = edges_attrib_dict.keys()
    Gcolors = edges_attrib_dict.values()

    #G.add_edges_from(edges_list) # this will only add edges but not the color info. which I might want later
    for (u,v) in edges_list:
        G.add_edge(u , v, color=edges_attrib_dict[(u,v)])

    # node attributes
    node_module_dict = readDict(node_attrib_file, 1, 7, '\t')
    for node in node_module_dict:
        G.add_node(node, moduleNumb = node_module_dict[node])

    # calculate avg. relative abundance of the node
    node_abundance = get_node_abundance(node_abund_file, delim)
    #and make a list and submit as node_size to draw()
    node_sizes_dict = node_weights_to_sizes(node_abundance)
    Gnode_sizes_list = [node_sizes_dict[i] for i in G.nodes()]

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
            color_ = 'cyan'
        samModG.add_edge(u , v, color=color_)
        #samModG_edges_list.append()
        samModGcolors.append(color_)
   
    #print samModGcolors

    # identify the category the OTU belongs to
    node_categ_dict = identify_node_categ(node_module_dict)
    # create a color map based on phyla
    categ_color_map = create_color_map(node_categ_dict.values())
    for node in samModG.nodes():
        samModG.add_node(node, category = node_categ_dict[node])

    samModGnode_sizes_list = [node_sizes_dict[i] for i in samModG.nodes()]
    #print samModGnode_sizes_list

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


    plt.title('OTUs colored as per modules. Intermodule edges allowed.')
    #nx.draw(G, edgelist=edges_list, edge_color = Gcolors, node_color = nodeColor, with_labels = False, node_size=Gnode_sizes_list)#, style='dashed')
    nx.draw(G, edgelist=edges_list, edge_color = Gcolors, node_color = nodeColor, with_labels = False)#, style='dashed')
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

    text_font = 'Arial' # 'Times' 'Helvetica'
    '''
    #http://stackoverflow.com/questions/18821795/how-can-i-get-list-of-font-familyor-name-of-font-in-matplotlib
    '''


    nx.draw(samModG, edgelist=edges_between_nodes_from_same_module, edge_color = samModGcolors, pos=pos, node_color = nodeColor, labels = new_labels, with_labels = True, node_size=samModGnode_sizes_list, font_size=8, font_family=text_font)
    #nx.draw(samModG, edgelist=edges_between_nodes_from_same_module, edge_color = samModGcolors, pos=pos, node_color = nodeColor, labels = new_labels, font_size=8, with_labels = True)
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

