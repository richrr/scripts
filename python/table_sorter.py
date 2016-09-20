import os
import sys
from utils import *
import operator
from time import localtime, strftime
import argparse
import math
import numpy
import matplotlib.pyplot as plt
import brewer2mpl
import pylab

#usage: python ~/scripts/python/table_sorter.py -i Group_otu_table_sorted_L2.txt

def main(args):

    parser = argparse.ArgumentParser(description='Sort table of taxonomic summary')
    parser.add_argument('-i', '--infile')
    parser.add_argument('-o', '--outfile', default="sortablefile.html")  # output filename
    parser.add_argument('-d', '--delimiter', default='\t')  # delimiter for file
    parser.add_argument('-f', '--float', action='store_true', default=False) # show numbers as decimal instead of percent (default)
    
    args = parser.parse_args()

    if len(sys.argv)==1 :
        parser.print_help()
        sys.exit('\natleast one argument required\n')

    infile  = args.infile
    outpfile = infile+args.outfile
    delim = args.delimiter
    showfloat = args.float
       
    lines = read_file(infile)
    if 'Taxon' not in lines[0]:
        sys.exit('\nheader line is required, tab-delimited, starts with taxon, followed by samples\n')

    content = ''
    content += FIXED_STRING()
    content += BEGIN_STRING(showfloat, infile)
    counter = 0
    TABS = '\t\t\t\t\t'
    for l in lines:
        counter += 1
        l = l.strip()
        if '#' in l or not l:
            continue
        if counter == 1:
            content += """\
            <table id="example" class="display" cellspacing="0" width="100%">
				<thead>
					<tr>

"""
	    val = [TABS+'<th>'+i+'</th>' for i in l.split(delim)]
            content = content + '\n'.join(val) + """\


					</tr>
				</thead>

"""
	    content += """\
				<tbody>

"""
	    continue
	content += """\
					<tr>
"""
        val = list()
        for i in l.split(delim):
            r = i
            try:
                if not showfloat:
                    i = float(i)*100
                r = "%.1f" % float(i)
            except:
                pass
            val.append(TABS+'<td>'+str(r)+'</td>')
        #val = [TABS+'<td>'+i+'</td>' for i in l.split(delim)]
	content = content + '\n'.join(val) + """\

					</tr>

"""
    content += END_STRING()
    #print content
    writeTXT_to_file(content, outpfile)


def END_STRING():
    string = """\
			
		</section>
	</div>

	
</body>
</html>

"""
    return string


def BEGIN_STRING(showfloat, infile):
    string = """\
<body class="dt-example">
	<div class="container">
		<section>
			<h1>Taxonomic summary<span> Sortable </span></h1>


"""
    string += '<h2>File: %s</h2>' %infile
    if not showfloat:
        string += '<h3>Numbers are in % relative abundance</h3>'
    return string

def FIXED_STRING():
    string = """\
<!DOCTYPE html>
<html>
<head>
	<meta charset="utf-8">
	<link rel="shortcut icon" type="image/ico" href="http://www.datatables.net/favicon.ico">
	<meta name="viewport" content="initial-scale=1.0, maximum-scale=2.0">

	<title>DataTables example - Default ordering (sorting)</title>
	<link rel="stylesheet" type="text/css" href="http://www.hort.vt.edu/microeco/RR/DataTables-1.10.2/media/css/jquery.dataTables.css">
	<link rel="stylesheet" type="text/css" href="http://www.hort.vt.edu/microeco/RR/DataTables-1.10.2/examples/resources/syntax/shCore.css">
	<link rel="stylesheet" type="text/css" href="http://www.hort.vt.edu/microeco/RR/DataTables-1.10.2/examples/resources/demo.css">
	<style type="text/css" class="init">

	</style>
	<script type="text/javascript" language="javascript" src="http://www.hort.vt.edu/microeco/RR/DataTables-1.10.2/media/js/jquery.js"></script>
	<script type="text/javascript" language="javascript" src="http://www.hort.vt.edu/microeco/RR/DataTables-1.10.2/media/js/jquery.dataTables.js"></script>
	<script type="text/javascript" language="javascript" src="http://www.hort.vt.edu/microeco/RR/DataTables-1.10.2/examples/resources/syntax/shCore.js"></script>
	<script type="text/javascript" language="javascript" src="http://www.hort.vt.edu/microeco/RR/DataTables-1.10.2/examples/resources/demo.js"></script>
	<script type="text/javascript" language="javascript" class="init">

$(document).ready(function() {
	$('#example').dataTable( {
		"order": [[ 0, "asc" ]]
	} );
} );

	</script>
</head>
"""
    return string


if __name__=='__main__':
    datetime = strftime("%a, %d %b %Y %I:%M:%S %p", localtime())
    cmd = 'echo ' + datetime 
    os.system(cmd)
    main(sys.argv)

