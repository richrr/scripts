import os
import sys
from utils import *
import operator
from time import localtime, strftime
import argparse
import re
from Bio import SeqIO
from abifpy import Trace
import pylab
from Bio.Blast import NCBIXML
from Bio.Blast import NCBIWWW
from Bio import Entrez

# usage: python /home/williamslab/scripts/python/blast_search_and_parse.py -i f.fasta

def main(args):
    # read sanger seq files
    # output to format
    # blast
    parser = argparse.ArgumentParser(description='Analyze Sanger seq abi files')
    parser.add_argument('-i', '--infile') # file containing blast results to be processed
    parser.add_argument('-o', '--outfile', default="default_output_file")  # output filename
    parser.add_argument('-d', '--delimiter', default='\t')  # delimiter for file
    parser.add_argument('-e', '--extension', default='fasta')  # output file format (extension) options allowed: fasta
    #parser.add_argument('-l', '--minseqlength', default=250, type=int) # minimum sequence length before trimming
    #parser.add_argument('-q', '--phredqthreshold', default=25, type=int) # phred quality threshold, 1 error in approx. 316 bases
    # http://www.sourcecodebrowser.com/python-biopython/1.59/namespace_bio_1_1_seq_i_o_1_1_abi_i_o.html#a7fb9b10aa5e88c3f1da16b212f364a3a
    #parser.add_argument('-r', '--revprimer', default='1492R')  # name of reverse primer
    #parser.add_argument('-s', '--fillerstring', default='SWG__') # string to make the header name >10 chars. long to make libshuff compatible
    
    args = parser.parse_args()

    if len(sys.argv)==1 :
        parser.print_help()
        sys.exit('\natleast one argument required\n')

    infile  = args.infile
    outfilename = args.outfile
    delim = args.delimiter
    extens = args.extension
    #segment = args.minseqlength


    #query_blast(infile, extens, outfilename, delim)
    query_blast_individ(infile, extens, outfilename, delim)

#http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc84

# blast all queries at once as a single string
def query_blast(infile, extens, outfilename, delim):
    xml_files = list()  # list of xml filenames of all queries # since only one query of multiple subqueries is sent, this will be only one file
    fasta_str = ''.join([rec.format(extens) for rec in SeqIO.parse(infile, extens)]) 
    #fasta_str='tacttgttgatattggatcgaacaaactggagaaccaacatgctcacgtcacttttagtcccttacatattcctc'
    results_handle = NCBIWWW.qblast("blastn", "nt", fasta_str ,\
            hitlist_size=10, nucl_penalty=-2, nucl_reward=1, megablast="TRUE")
    xml_file = outfilename + "_megablast_results.xml"
    xml_files.append(xml_file) # add file name to list of xmls
    save_file = open(xml_file, "w")
    save_file.write(results_handle.read())
    save_file.close()
    results_handle.close()
    parse_blast_xmls(outfilename, xml_files, delim)

# blast each query by itself (individually)
def query_blast_individ(infile, extens, outfilename, delim):
    xml_files = list()  # list of xml filenames of all queries 
    for rec in SeqIO.parse(infile, extens):
        fasta_str = rec.format(extens)
        results_handle = NCBIWWW.qblast("blastn", "nt", fasta_str ,\
            hitlist_size=10, nucl_penalty=-2, nucl_reward=1, megablast="TRUE")
        xml_file = outfilename + "_query_" + rec.id + "_individual_megablast_results.xml"
        xml_files.append(xml_file) # add file name to list of xmls
        save_file = open(xml_file, "w")
        save_file.write(results_handle.read())
        save_file.close()
        results_handle.close()
    parse_blast_xmls(outfilename, xml_files, delim)


'''
qblast options:
 http://www.ncbi.nlm.nih.gov/BLAST/Doc/urlapi.html

why to use megablast: 
 http://www.ncbi.nlm.nih.gov/blast/Why.shtml 
 http://www.ncbi.nlm.nih.gov/staff/tao/URLAPI/new/node81.html

http://www.ncbi.nlm.nih.gov/books/NBK25497/#chapter2.T._entrez_unique_identifiers_ui
http://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.EFetch

Other non-python resources if needed:
https://www.biostars.org/p/13452/
http://www.polarmicrobes.org/?p=759
http://search.cpan.org/~cjfields/BioPerl/Bio/DB/Taxonomy.pm
http://search.cpan.org/~motif/Bio-LITE-Taxonomy-NCBI-0.09/lib/Bio/LITE/Taxonomy/NCBI.pm
'''

#Usage, opens an outfile and then parses any number of .xml files (from the xml_files list) into that outfile, printing all hits
def parse_blast_xmls(outfilename, xml_files, delim):
    outfile = outfilename+"_parsed_megablast_results.txt"
    OUT = open(outfile, 'w')
    ########## Double check, THESE RESULTS MAY NOT SORTED FROM TOP TO WORST HIT #################
    OUT.write("Query Name\tQuery Length\tSubject Name\tSubject Length\tAlignment Length\tQuery Sequence\tHsp Score\tHsp Expect\tHsp Identities\tPercent Match\tNumber_of_gaps\n")
    for xml_file in xml_files:
        result_handle = open(xml_file)
        blast_records = NCBIXML.parse(result_handle)
        for rec in blast_records:
           for alignment in rec.alignments:
              for hsp in alignment.hsps:
                 OUT.write(str(rec.query) + '\t' + str(rec.query_length) + '\t' + str(alignment.title) + '\t' + str(alignment.length) + '\t' + str(hsp.align_length) + '\t' + str(hsp.query) + '\t' + str(hsp.score) + '\t' + str(hsp.expect) + '\t' + str(hsp.identities) + '\t' + str(float(hsp.identities)/int(hsp.align_length)) + '\t' + str(hsp.gaps) + '\n')
        result_handle.close()
    OUT.close()
    parse_gis_megab_result(outfilename, outfile, delim)


# http://www.ncbi.nlm.nih.gov/books/NBK25497/table/chapter2.T._entrez_unique_identifiers_ui/?report=objectonly
# https://www.biostars.org/p/99855/

def parse_gis_megab_result(outfilename, parsed_blast_infile, delim):
    # for each query seq. in the file, use the gi for top three (assuming the provided results were already desc. sorted as per top hits)
    lines = read_file(parsed_blast_infile)
    queries_gis = dict()
    gid_col = 0  # 0 based index
    perc_match_col = 0 # 0 based index
    for l in lines:
        l = l.strip()
        if '#' in l or not l:
            continue
        contents = l.split(delim)
        if 'Subject Name' in contents or 'Percent Match' in contents:
            gid_col =  contents.index('Subject Name')
            perc_match_col = contents.index('Percent Match')
            continue
        query = contents[0]
        hit_id = contents[gid_col].split('|')[1] # gi|703126872|ref|XM_010105390.1| Morus notabilis hypothetical protein partial mRNA
        perc_match = contents[perc_match_col]
        if query not in queries_gis:
            queries_gis[query] = [[hit_id, perc_match]] # list of lists
            #queries_gis[query] = [hit_id]
        elif len(queries_gis[query]) < 3:
            #queries_gis[query].append(hit_id)
            queries_gis[query].append([hit_id, perc_match])
    #print queries_gis
    #download_parse_gis_per_query(outfilename, queries_gis)
    download_parse_indiv_gis(outfilename, queries_gis)

def download_parse_indiv_gis(outfilename, queries_gis):
    results_organism_taxa_list = list()
    header_Str = '#query\tgi\tperc_match\tGBSeq_primary-accession\tGBSeq_organism\tGBSeq_taxonomy'
    for q in queries_gis:
      #for gid in queries_gis[q]:
      for gid_pm in queries_gis[q]:
        gid = gid_pm[0]
        pm = gid_pm[1]
        
        filename = "gb_" + gid +".xml"
        Entrez.email = "richrr@vt.edu"
        
        # avoid redownloading same files
        if not os.path.isfile(filename):
            handle = Entrez.efetch(db="nucleotide", id=gid, rettype="gb", retmode="xml")

            #save entrez result to xml file
            out_handle = open(filename, "w")
            out_handle.write(handle.read())
            out_handle.close()
            handle.close()

        # parse the entrez xml file
        gb_records = Entrez.parse(open(filename), validate=False)
        for record in gb_records:
            #str = q + '\t' + gid + '\t' + record['GBSeq_primary-accession'] + '\t' + record['GBSeq_organism'] + '\t' + record['GBSeq_taxonomy']
            str = q + '\t' + gid + '\t' + pm + '\t' + record['GBSeq_primary-accession'] + '\t' + record['GBSeq_organism'] + '\t' + record['GBSeq_taxonomy']
            results_organism_taxa_list.append(str)
            #print "%s" % (str)

    results_taxa_file = outfilename+"_parsed_wtaxa_results.txt"
    writeLIST_to_file(results_organism_taxa_list, results_taxa_file)


def download_parse_gis_per_query(outfilename, queries_gis):
    # this would reduce the number of calls to NCBI by simultaneously downloading file for multiple gis
    # need to check how it behaves if the same gid is given multiple times in a list for a query
    # may not be able to print the percent match associated with each gi for a query
    # a hack would be to search the gi name from the record, and use the corressp perc match from the list
    # incase the gi is repeated in the list, then the perc match really doesn't matter
    results_organism_taxa_list = list()
    for q in queries_gis:
        #gid_list = queries_gis[q]
        gid_list = [gid_pm[0] for gid_pm in queries_gis[q]]
        filename = "q_" + q +".xml"
        Entrez.email = "richrr@vt.edu"
        handle = Entrez.efetch(db="nucleotide", id=gid_list, rettype="gb", retmode="xml")

        #save entrez result to xml file
        out_handle = open(filename, "w")
        out_handle.write(handle.read())
        out_handle.close()
        handle.close()

        # parse the entrez xml file
        gb_records = Entrez.parse(open(filename), validate=False)
        for record in gb_records:
            str = q + '\t' + record['GBSeq_other-seqids'][1] + '\t' + record['GBSeq_primary-accession'] + '\t' + record['GBSeq_organism'] + '\t' + record['GBSeq_taxonomy']
            results_organism_taxa_list.append(str)
            #print "%s" % (str)

    results_taxa_file = outfilename+"_parsed_gis_per_query_wtaxa_results.txt"
    writeLIST_to_file(results_organism_taxa_list, results_taxa_file)
    


'''
gid="186972394"
gid="568824607"
filename = "gb_" + gid+".xml"
Entrez.email = "richrr@vt.edu"
handle = Entrez.efetch(db="nucleotide", id=gid, rettype="gb", retmode="xml")

#save to xml file
out_handle = open(filename, "w")
out_handle.write(handle.read())
out_handle.close()
handle.close()

# parse the entrez xml file
gb_records = Entrez.parse(open(filename), validate=False)
for record in gb_records:
    print "%s\t%s\t%s" % (record['GBSeq_primary-accession'], record['GBSeq_organism'] , record['GBSeq_taxonomy'])

'''
'''
http://biopython.org/DIST/docs/tutorial/Tutorial.html#sec:entrez-webenv
use the taxonomic strategy if the above gi identifiers fail
>>> from Bio import Entrez
>>> Entrez.email = "A.N.Other@example.com"     # Always tell NCBI who you are
>>> handle = Entrez.esearch(db="Taxonomy", term="Cypripedioideae")
>>> record = Entrez.read(handle)
>>> record["IdList"][0]
'158330'

Now, we use efetch to download this entry in the Taxonomy database, and then parse it:

>>> handle = Entrez.efetch(db="Taxonomy", id="158330", retmode="xml")
>>> records = Entrez.read(handle)

Again, this record stores lots of information:

>>> records[0].keys()
[u'Lineage', u'Division', u'ParentTaxId', u'PubDate', u'LineageEx',
 u'CreateDate', u'TaxId', u'Rank', u'GeneticCode', u'ScientificName',
 u'MitoGeneticCode', u'UpdateDate']

We can get the lineage directly from this record:

>>> records[0]["Lineage"]
'cellular organisms; Eukaryota; Viridiplantae; Streptophyta; Streptophytina;
 Embryophyta; Tracheophyta; Euphyllophyta; Spermatophyta; Magnoliophyta;
 Liliopsida; Asparagales; Orchidaceae'

'''

if __name__=='__main__':
    datetime = strftime("%a, %d %b %Y %I:%M:%S %p", localtime())
    cmd = 'echo ' + datetime 
    os.system(cmd)
    main(sys.argv)



