import os
import sys
from utils import *
import re
import itertools
import numpy as np
from collections import Counter
import operator


### create a dict of original to comma replaced seeds entries. Would make life much easier.

'''
###### usage #####
cd /nfs3/PHARM/Morgun_Lab/richrr/Type2_Diabetes/microbes-regulating-mmp12/OV-selection-sep-2019/workspace/m-m/merged/corr/st.t2.OTU_176957_seed_seed/p1_cp1_cfdr1.1/per_analysis

python /nfs3/PHARM/Morgun_Lab/richrr/scripts/python/parse-seed-for-genome-id-updated.py /nfs3/PHARM/Morgun_Lab/richrr/Type2_Diabetes/microbes-regulating-mmp12/OV-selection-sep-2019/workspace/m-m/merged/corr/st.t2.16S_shotgun_seed/p1_cp1_cfdr1.1/per_analysis/OTU_176957.seed.0.70.thresh.csv ~/Morgun_Lab/richrr/db/SEED/seed.fa ~/Morgun_Lab/richrr/db/SEED/seed-ids-only.fa



# testing only ##
cd /nfs3/PHARM/Morgun_Lab/richrr/Type2_Diabetes/microbes-regulating-mmp12/OV-selection-sep-2019/workspace/m-m/merged/corr/st.t2.OTU_176957_seed_seed/

grep "2'.3'-cyclic-nucleotide 2'-phosphodiesterase (EC 3.1.4.16)" ~/Morgun_Lab/richrr/Type2_Diabetes/DNAShotgunMetagenome/analysis/t2-corr/merged/otu_seed/p1_cp0.05_cfdr0.1/per_analysis/all_samples/ip_0.3_cp_0.1_comb-fdr_1/lactis-seed-corr-0.6thresh.csv > del.txt

grep "2-succinyl-6-hydroxy-2.4-cyclohexadiene-1-carboxylate synthase (EC 4.2.99.20)" ~/Morgun_Lab/richrr/Type2_Diabetes/DNAShotgunMetagenome/analysis/t2-corr/merged/otu_seed/p1_cp0.05_cfdr0.1/per_analysis/all_samples/ip_0.3_cp_0.1_comb-fdr_1/lactis-seed-corr-0.6thresh.csv >> del.txt


python /nfs3/PHARM/Morgun_Lab/richrr/scripts/python/parse-seed-for-genome-id-updated.py del.txt ~/Morgun_Lab/richrr/db/SEED/seed.fa ~/Morgun_Lab/richrr/db/SEED/seed-ids-only.fa

'''



infile= sys.argv[1]
refseedfile = sys.argv[2]

refseedNAMES = sys.argv[3]


n = 80
seed_entries = list()
seed = read_file(infile)

#'''
# use this if creating the seed entries by splitting the seed levels
# some of these may not be sanitized, i.e. they may have '.' instead of ','
# and may not end up being detected
delim = ';'

for line in seed:
    line=line.strip()
    if not line:
        continue
    #print line
    conts = line.split(delim)
    seed_entries.append(conts[-2])
    #sys.exit(0)

writeLIST_to_file(seed_entries, "seed-entries-unsanitized.txt")
#writeLIST_to_file(seed_entries, "seed-entries.txt")
#'''

# use this if directly using pre-sanitized seed entries
#seed_entries = [line.strip() for line in seed]
#print seed_entries




# get the list of ids from ref file
'''
refSeedNamesLIST = list() 
with open(refseedfile) as f:
    for line in f:
        line=line.strip()
        if line.startswith(">"):
             refSeedNamesLIST.append(line)  
'''
refSeedNamesLIST = read_file(refseedNAMES)


# this section only tries to santize the problematic seed entries to match the ones in the reference
unsantized_cant_find_seeds = list()
santized_can_find_seeds = list()

n = 0  # counter of how many seed terms processed
for substring in seed_entries:
  n = n+1
  print n
  if any(substring in line for line in refSeedNamesLIST):
      santized_can_find_seeds.append(substring)
      continue
  else:
      print ("Cannot find: " + substring + "   Trying something.")
      
      str = substring
      
      # try if a simple replace of dot fixes the issue, no worry about (...)
      finalstr = re.sub(r'\.', ',' , str)
      if any(finalstr in line for line in refSeedNamesLIST):
          santized_can_find_seeds.append(finalstr)
          print ("Found it in try 1")
          continue
      
      
      if '(' in str and ')' in str: 
        # get the keep as is pattern inside ()
        patt = re.findall(r"\(.+\)", str)[-1]  #if there are many patterns with ()
        # replace . with , in str
        newstr = re.sub(r'\.', ',' , str)
        # replace the (....) with keep as is pattern in the newstr
        finalstr = re.sub(r"\(.+\)", patt , newstr)
        #print finalstr

        if any(finalstr in line for line in refSeedNamesLIST):
          santized_can_find_seeds.append(finalstr)
          print ("Found it in try 2")
          continue

      print substring
      unsantized_cant_find_seeds.append(substring)
      

if len(seed_entries) == len(santized_can_find_seeds):
    print "\nFound ALL"
else:
    print "Couldn't find these: Ignoring them and moving on"
    print(unsantized_cant_find_seeds)
    print(len(unsantized_cant_find_seeds))
    #sys.exit(0)



'''
import re
str="2'.3'-cyclic-nucleotide 2'-phosphodiesterase (EC 3.1.4.16)"
# get the keep as is pattern inside ()
patt = re.findall(r"\(.+\)", str)[0]

# replace . with ,
newstr = re.sub(r'\.', ',' , str)

# go back and replace the (....) with keep as is pattern
finalstr = re.sub(r"\(.+\)", patt , newstr)
finalstr
'''



#counter = 0

outfile = list()
seq=''
with open(refseedfile) as f:
    for line in f:
        #counter = counter + 1
        #if counter == 10000:
        #    break
        
        line=line.strip()
        #print line
        
        if line.startswith(">"):
            if any(substring in line for substring in santized_can_find_seeds): #https://stackoverflow.com/questions/8122079/python-how-to-check-a-string-for-substrings-from-a-list
                outfile.append(line)
                #print(line)
                flag = 1
            else:
                flag = 0
            
            strs = [seq[i:i+n] for i in range(0, len(seq), n)]
            # do not add empty lines
            if len(strs) > 0:
                outfile.extend(strs)
            
            seq=''
            continue


        if flag == 1:
            seq = seq+line
        
writeLIST_to_file(outfile, "seed-seqs.txt")




# parse the results to count genomes

seed_seq_ids = list() # file containing only the full seed header, no seqs
seed_genome_dict = dict() # file containing a dictionary of the seed term for which the hit was found

with open("seed-seqs.txt") as f:
    for line in f:
        line=line.strip()
        if line.startswith(">"):
            seed_seq_ids.append(line)
            
            for substring in santized_can_find_seeds:
                if substring in line:
                    #print substring
                    # add to dict
                    if substring in seed_genome_dict:
                        val = seed_genome_dict[substring]
                        val.append(line)
                        seed_genome_dict[substring] = val
                    else:
                        seed_genome_dict[substring] = [line]
                    # not breaking loop to allow multiple matches to protein

writeLIST_to_file(seed_seq_ids, "seed-seqs-ids.txt")
writeDICT_to_file(seed_genome_dict, "seed-genome-dict.txt")




# extract the taxa name from the []
def parse_genomes(gs):
    genomes = list()
    for g in gs:
        #print g
        matches = re.findall(r"\[(.+)\]", g)
        for match in matches:
            #print match
            genomes.append(match)
    return genomes


# this contains the actual taxa of the bacteria containing the seed protein
seed_taxa_dict = dict() # note that the same genome can have many occurrences of the seed
seed_utaxa_dict = dict() # unique taxa, so we can later count which taxa contains how many of our seeds
for k,v in seed_genome_dict.items():
    newv = parse_genomes(v)
    seed_taxa_dict[k] = newv
    seed_utaxa_dict[k] = list(set(newv))

writeDICT_to_file(seed_taxa_dict, "seed-taxa-dict.txt")
writeDICT_to_file(seed_utaxa_dict, "seed-uniq-taxa-dict.txt")

#print seed_utaxa_dict.keys()
writeLIST_to_file(seed_utaxa_dict.keys(), "seed-dict-keys.txt")






# count which taxa contains how many seeds:


def sort_taxa_occurence(somelistoflists, outfile):
    # merges list of list into one list # https://stackoverflow.com/questions/716477/join-list-of-lists-in-python
    one_list = list(itertools.chain.from_iterable(somelistoflists)) 
    # https://stackoverflow.com/questions/46580499/count-occurrences-of-an-element-in-list
    taxa_occurences =  Counter(one_list)
    #writeDICT_to_file(taxa_occurences , outfile)

    # reverse sorted as per occurrences # https://stackoverflow.com/questions/613183/how-to-sort-a-dictionary-by-value
    sorted_x = sorted(taxa_occurences.items(), key=operator.itemgetter(1), reverse=True)
    #print sorted_x

    #https://stackoverflow.com/questions/3820312/python-write-a-list-of-tuples-to-a-file
    with open(outfile, 'w') as fp:
        fp.write('\n'.join('%s\t%s' % x for x in sorted_x))



# each taxa occurence indicates one seed. so the score indicates number of seeds covered by taxa
sort_taxa_occurence(seed_utaxa_dict.values(), 'uniq-taxa-covering-seeds.txt')


# taxa can occur multiple times per seed. so the score can be greater than number of seeds entries for this code
sort_taxa_occurence(seed_taxa_dict.values(), 'taxa-covering-seeds.txt')


print "Done"



