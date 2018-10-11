import sys
import os
import networkx as nx
from utils import *
import pandas as pd


#usage: python ~/Morgun_Lab/richrr/scripts/python/get_microbiota_regulated_neighbors.py /nfs3/PHARM/Morgun_Lab/richrr/Type2_Diabetes/RNA-Seq/analysis/t2-multi-organ/gene_pairs.txt ~/Morgun_Lab/richrr/Type2_Diabetes/RNA-Seq/t2-rna-seq-cyto/multi-organ-nets/spf-gf-diet-degs-microbe-diet-interaction.txt

infile = sys.argv[1]
spf_gf_diet_attrib_file = sys.argv[2]

#print infile
#print spf_gf_diet_attrib_file


df = pd.read_csv(spf_gf_diet_attrib_file, sep='\t', index_col=0)
#print(df)

concorfc_fdr = 'Concordant FC and interaction FDR < 15%'
discorfc_fdr = 'Discordant FC and interaction FDR < 15%'

#print(df['Symbol']['ENSMUSG00000030643_I'])
#print(df['Discordant FC and interaction FDR < 15%']['ENSMUSG00000030643_I'])



G=nx.read_edgelist(infile, nodetype=str, delimiter='\t')

pairs = list(G.edges)
#print(pairs)
print(len(pairs))

genes = list(G.nodes)
#print(genes)
print(len(genes))

outputlist = ['\t'.join(['Node', 'Numb_Neighbor_Genes', 'Neighbor_Genes_w_Discordant_HfvsNcdFC_due_to_Microbes_and_interactionFDR<15%' , 'Neighbor_Genes_w_Concordant_HfvsNcdFC_due_to_Microbes_and_interactionFDR<15%' ,'Neighbor_Genes'])]

for n in genes:
    ngbrs = list(G.adj[n])
    #print n , G.degree[n], ngbrs
    dissigcounter = 0
    consigcounter = 0
    for g in ngbrs:
        discordant = df[discorfc_fdr][g]
        concordant = df[concorfc_fdr][g]
        
        dissigcounter = dissigcounter + discordant
        consigcounter = consigcounter + concordant
        
        #print g, discordant, concordant
    print n , G.degree[n], dissigcounter, consigcounter, ','.join(ngbrs)
    outstr = '\t'.join([n , str(G.degree[n]), str(dissigcounter), str(consigcounter), ','.join(ngbrs)])
    outputlist.append(outstr)


writeLIST_to_file(outputlist, 'microbiota_regulated_neighbors_of_genes_output.txt')


