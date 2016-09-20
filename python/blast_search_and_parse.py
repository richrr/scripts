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
from xml.dom import minidom

# usage: python /home/williamslab/scripts/python/blast_search_and_parse.py -i f.fasta

#http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc84

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
    parser.add_argument('-b', '--blastoption', default='megablast')  # by default, run megablast. Other option allowed is blastn
    #parser.add_argument('-s', '--fillerstring', default='SWG__') # string to make the header name >10 chars. long to make libshuff compatible
    
    args = parser.parse_args()

    if len(sys.argv)==1 :
        parser.print_help()
        sys.exit('\natleast one argument required\n')

    infile  = args.infile
    outfilename = args.outfile
    delim = args.delimiter
    extens = args.extension
    blast_option = args.blastoption


    query_blast(infile, extens, outfilename, delim, blast_option)
    #parse_gis_megab_result(outfilename, "se_v0_july_fwd_dir_qtrim_phred_25_parsed_megablast_results.txt", delim, blast_option)



def parse_blast_to_check_error(xml_file, blast_option) :
    blast_error_list = ["Query 'lcl|60945 stage_v0_RC_1KXD137RV0+1492R_reverse_complement' (# 4): Warning: Could not calculate ungapped Karlin-Altschul parameters due to an invalid query sequence or its translation. Please verify the query sequence(s) and/or filtering options" ,\
            "[blastsrv4.REAL]: Error: CPU usage limit was exceeded, resulting in SIGXCPU (24)."]
    rerun = "False"
    xmldoc = minidom.parse(xml_file)
    itemlist = xmldoc.getElementsByTagName('Iteration_message') 
    for s in itemlist :
        err_msg = s.childNodes[0].nodeValue
        if err_msg in blast_error_list:
            print "Need to rerun %s individually on seqs. in this query file because\n %s" % (blast_option, err_msg)
            rerun = "True"
            break
    return rerun


def batch_iterator(iterator, batch_size) :
    """Returns lists of length batch_size.
 
    This can be used on any iterator, for example to batch up
    SeqRecord objects from Bio.SeqIO.parse(...), or to batch
    Alignment objects from Bio.AlignIO.parse(...), or simply
    lines from a file handle.
 
    This is a generator function, and it returns lists of the
    entries from the supplied iterator.  Each list will have
    batch_size entries, although the final list may be shorter.
    """
    entry = True #Make sure we loop once
    while entry :
        batch = []
        while len(batch) < batch_size :
            try :
                entry = iterator.next()
            except StopIteration :
                entry = None
            if entry is None :
                #End of file
                break
            batch.append(entry)
        if batch :
            yield batch



# submit seqs. in batches, for each batch: blast all queries at once as a single string
def query_blast(infile, extens, outfilename, delim, blast_option):
    xml_files = list()  # list of xml filenames of all queries # since only one query of multiple subqueries is sent, this will be only one file
    record_iter = SeqIO.parse(open(infile), extens)
    for i, batch in enumerate(batch_iterator(record_iter, 10)) :
        filename = "%s_group_%i.%s" % (outfilename, (i+1), extens)
        handle = open(filename, "w")
        count = SeqIO.write(batch, handle, extens)
        handle.close()
        print "Wrote %i records to %s" % (count, filename)

        fasta_str = ''.join([rec.format(extens) for rec in SeqIO.parse(filename, extens)])
        results_handle = None
        if blast_option == 'megablast':
            results_handle = NCBIWWW.qblast("blastn", "nt", fasta_str ,\
                   expect=10, hitlist_size=10, nucl_penalty=-2, nucl_reward=1, megablast="TRUE")
        elif blast_option == 'blastn':
            results_handle = NCBIWWW.qblast("blastn", "nt", fasta_str ,\
                   expect=10, hitlist_size=10) # nucl_penalty=-3, nucl_reward=2, megablast="FALSE", gapcosts = '5 2' # gap open and extend cost
        xml_file = filename + "_" + blast_option + "_results.xml"
        save_file = open(xml_file, "w")
        save_file.write(results_handle.read())
        save_file.close()
        results_handle.close()
        # parse the xml_file for errors
        rerun = parse_blast_to_check_error(xml_file, blast_option)
        # if returns error (rerun True), rerun the seqs in that chunk (filename) as separate queries
        # not checking the individual runs for errors
        if rerun == 'True':
            print "Rerunning seqs. from %s as individual queries" %filename
            indivi_query_xml_files = query_blast_individ(filename, extens, outfilename, delim, blast_option)
            xml_files.extend(indivi_query_xml_files)
        else:
            xml_files.append(xml_file) # add file name to list of xmls

    # at this point, hopefully all queries have been successfully blasted
    parse_blast_xmls(outfilename, xml_files, delim, blast_option)    


# blast each query by itself (individually)
def query_blast_individ(infile, extens, outfilename, delim, blast_option):
    xml_files = list()  # list of xml filenames of all queries 
    for rec in SeqIO.parse(infile, extens):
        fasta_str = rec.format(extens)
        results_handle = None
        if blast_option == 'megablast':
            results_handle = NCBIWWW.qblast("blastn", "nt", fasta_str ,\
                  hitlist_size=10, nucl_penalty=-2, nucl_reward=1, megablast="TRUE")
        elif blast_option == 'blastn':
            results_handle = NCBIWWW.qblast("blastn", "nt", fasta_str ,\
                   expect=10, hitlist_size=10) # nucl_penalty=-3, nucl_reward=2, megablast="FALSE", gapcosts = '5 2' # gap open and extend cost
        xml_file = infile + "_query_" + rec.id + "_indiv_" + blast_option + "_results.xml"
        xml_files.append(xml_file) # add file name to list of xmls
        save_file = open(xml_file, "w")
        save_file.write(results_handle.read())
        save_file.close()
        results_handle.close()
    return xml_files



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
def parse_blast_xmls(outfilename, xml_files, delim, blast_option):
    outfile = outfilename+"_parsed_" + blast_option + "_results.txt"
    OUT = open(outfile, 'w')
    ########## Double check, THESE RESULTS MAY NOT SORTED FROM TOP TO WORST HIT #################
    OUT.write("Hit number\tQuery Name\tQuery Length\tHSP Id\tHSP Name\tHSP Title\tHSP Length\tAlignment Length\tQuery Sequence\tHsp Score\tHsp Expect\tHsp Identities\tPercent Match\tNumber_of_gaps\n")
    for xml_file in xml_files:
        result_handle = open(xml_file)
        blast_records = NCBIXML.parse(result_handle)
        for rec in blast_records:
           #for alignment in rec.alignments:
           for hit_num, alignment in enumerate(rec.alignments):
              for hsp in alignment.hsps:
                 OUT.write(str((hit_num+1)) + '\t' + str(rec.query) + '\t' + str(rec.query_length) + '\t' + \
                   str(alignment.hit_id) + '\t' + str(alignment.hit_def) + '\t' + str(alignment.title) + '\t' + \
                   str(alignment.length) + '\t' + str(hsp.align_length) + '\t' + str(hsp.query) + '\t' + \
                   str(hsp.score) + '\t' + str(hsp.expect) + '\t' + str(hsp.identities) + '\t' + \
                   str(float(hsp.identities)/int(hsp.align_length)) + '\t' + str(hsp.gaps) + '\n')
        result_handle.close()
    OUT.close()
    parse_gis_megab_result(outfilename, outfile, delim, blast_option)

# alignment.hit_id, alignment.hit_def, alignment.hit_accession
#  alignment.hit_num does not work

# http://www.ncbi.nlm.nih.gov/books/NBK25497/table/chapter2.T._entrez_unique_identifiers_ui/?report=objectonly
# https://www.biostars.org/p/99855/

def parse_gis_megab_result(outfilename, parsed_blast_infile, delim, blast_option):
    # for each query seq. in the file, use the gi for top three (assuming the provided results were already desc. sorted as per top hits)
    lines = read_file(parsed_blast_infile)
    queries_gis = dict()
    gid_col = 0  # 0 based index
    perc_match_col = 0 # 0 based index
    hitnum_col = 0
    query_name_col = 0
    for l in lines:
        l = l.strip()
        if '#' in l or not l:
            continue
        contents = l.split(delim)
        if 'HSP Id' in contents or 'Percent Match' in contents:
            hitnum_col = contents.index('Hit number')
            query_name_col = contents.index('Query Name')
            gid_col =  contents.index('HSP Id')
            perc_match_col = contents.index('Percent Match')
            continue
        query = contents[query_name_col]
        hit_id = contents[gid_col].split('|')[1] # gi|703126872|ref|XM_010105390.1
        perc_match = contents[perc_match_col]
        hit_num = contents[hitnum_col]
        if query not in queries_gis:
            queries_gis[query] = [[hit_id, perc_match, hit_num]] # list of lists
            #queries_gis[query] = [hit_id]
        elif len(queries_gis[query]) < 3: # if a hit has multiple hsps, this would take top 3 hsp from the same hit, not top 3 hits
          curr_count = len(queries_gis[query]) - 1 # 0 based index
          if queries_gis[query][curr_count][2] < hit_num :
             queries_gis[query].append([hit_id, perc_match, hit_num])
    #print queries_gis
    #download_parse_gis_per_query(outfilename, queries_gis)
    download_parse_indiv_gis(outfilename, queries_gis, blast_option)

def download_parse_indiv_gis(outfilename, queries_gis, blast_option):
    results_organism_taxa_list = list()
    header_Str = '#Hit_num\tQuery\tgi\tperc_match\tGBSeq_primary-accession\tGBSeq_organism\tGBSeq_taxonomy'
    results_organism_taxa_list.append(header_Str)
    for q in queries_gis:
      #for gid in queries_gis[q]:
      for gid_pm in queries_gis[q]:
        gid = gid_pm[0]
        pm = gid_pm[1]
        hnum = gid_pm[2]
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
            str = hnum + '\t' + q + '\t' + gid + '\t' + pm + '\t' + record['GBSeq_primary-accession'] + '\t' + record['GBSeq_organism'] + '\t' + record['GBSeq_taxonomy']
            results_organism_taxa_list.append(str)
            #print "%s" % (str)

    results_taxa_file = outfilename + "_" + blast_option + "_parsed_wtaxa_results.txt"
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



