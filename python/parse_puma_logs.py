import os
import sys
from utils import *
import operator
from time import localtime, strftime
import argparse
import re
import sh



def parse_count_files(f, delim='\t'):
    res = list()
    cont = read_file(f)
    #res.append(str(len(cont)))  # unique rows; however this might also contain empty rows
    vals = list()
    #print cont
    for l in cont:
        l = l.strip('\n')
        if l: # accounts for empty line
            #print l
            c = l.split(delim)
            vals.append(int(c[1].strip()))
    res.append(str(len(vals)))  # unique non-empty rows
    #print vals
    total = sum(vals)
    print total
    res.append(str(total))  # sum of values
    return res

def main(args):
    # extract the statistics from the log
    parser = argparse.ArgumentParser(description='extract the statistics from the log of code for metagenomic analysis')
    parser.add_argument('-i', '--indir')
    parser.add_argument('-o', '--outfile', default="puma-statistics.txt")  # output filename
    parser.add_argument('-d', '--delimiter', default='\t')  # delimiter for file
    
    args = parser.parse_args()

    if len(sys.argv)==1 :
        parser.print_help()
        sys.exit('\natleast one argument required\n')

    indir  = args.indir
    outfile = args.outfile
    delim = args.delimiter

    contents = sh.ls(sh.glob(indir+"*fastq*dir/run_info*log"))
    contents = [words for segments in contents for words in segments.split()] # split the elements in contents by space and put all in list
    # the above is useful if the ls command lists multiple files on a single line which is treated as a single item in the contents list
    print contents
    #sys.exit()

	# the numbers from megan's search are not the final outputs. Depending on megan's options, the outputs may vary, so specifically use the output files
	# to get the statistics
    statistics_list = [["Sample", "Total reads processed", "Reads with adapters" , "Reads that were too short" , 
                        "Reads written (passing filters)", "Total basepairs processed" , "Total written (filtered)",
                        "Number of reads","Low complexity","With valid hits","With SEED-ids","With KEGG-ids","Number of taxa identified",
                        "Number of SEED classes identified","Number of KEGG classes identified","Data processor required","Total reads",
                        "Assigned reads","Unassigned reads","Reads with no hits","Reads low comp.","KEGG rows", "KEGG values", "SEED rows", "SEED values", "TAXON rows", "TAXON values"]]

    for log in contents:
        #print "\nFile: ", log , "\n----After cutadapt----\n"
        print log
        log = log.strip()
        grep_cutadapt = sh.grep.bake("-A 9")
        cutadapt_ = grep_cutadapt("Summary", log)
        #print cutadapt_
        
        local_list = [log]
        for c in cutadapt_:
            c = c.strip()
            if len(c.split(':')) > 1: 
                numbs = [r.lstrip() for r in c.split(':')]
                local_list.append(numbs[1])

        
        #print "----After MEGAN----\n"
        grep_megan = sh.grep.bake("-A 16")

        megan_ = list()
        try:
            megan_ = grep_megan("Analyzing all matches", log)
            #megan_ = grep_megan("\'Analyzing all matches\'", log)
        except Exception, e:
            print "Found nothing for megan in log"
            print e.stderr
            
        
        #print megan_
        
        if len(megan_) > 0:
          for c in megan_:
            c = c.strip()
            if len(c.split(':')) > 1: 
                numbs = [r.lstrip() for r in c.split(':')]
                local_list.append(numbs[1])

        
        # count the number of rows and sum values for the kegg, seed, and taxon files.
        kegg_file = re.sub(r'/run_info\.\d+\.log', "/keggpath_count.txt", log)
        print kegg_file
        local_list.extend(parse_count_files(kegg_file))

        seed_file = re.sub(r'/run_info\.\d+\.log', "/seedpath_count.txt", log)
        print seed_file
        local_list.extend(parse_count_files(seed_file))

        taxon_file = re.sub(r'/run_info\.\d+\.log', "/taxonpath_count.txt", log)
        print taxon_file
        local_list.extend(parse_count_files(taxon_file))

        statistics_list.append(local_list)


    pr_str = ''
    for line in statistics_list:
        pr_str = pr_str + '\n' + delim.join(line)

    #print pr_str
    writeTXT_to_file(pr_str, outfile)

if __name__=='__main__':
    datetime = strftime("%a, %d %b %Y %I:%M:%S %p", localtime())
    cmd = 'echo ' + datetime 
    os.system(cmd)
    main(sys.argv)

