import urllib2
import sys
import os
from utils import *
import glob
import xml.etree.ElementTree as ET

to_do = sys.argv[2]

def get_files(patterns):
    for p in patterns:
        p = p.strip()
        if not p:
            continue
        outfile = p + ".xml"
        print outfile
        url = "https://www.ebi.ac.uk/ena/data/view/" + p + "&display=xml&download=xml&filename=" + outfile
        s = urllib2.urlopen(url)
        contents = s.read()
        file = open(outfile, 'w')
        file.write(contents)
        file.close()

def parse_files(files):
    resul = ['\t'.join(["File", "SecAccession", "Accession", "DKDC1", "ExpAccession", "SampleID", "DKDC2", "Age", "HPV_Status", "HPV_Status_At_Followup", "Tissue_Type", "ENA_Checklist", "ENA_Spot_Count", "ENA_Base_Count"])]
    for f in files:
        print f
        if not f:
            continue
        local = [f]
        tree = ET.parse(f)
        root = tree.getroot()
        #for ids in root.findall('SAMPLE'):
        #    print ids.tag, ids.attrib
        prim_id = root.find('./SAMPLE/IDENTIFIERS/PRIMARY_ID').text
        ext_id = root.find('./SAMPLE/IDENTIFIERS/EXTERNAL_ID').text
        #print prim_id, ext_id
        local.extend([prim_id, ext_id])
        
        sample_info = root.findall('./SAMPLE/SAMPLE_LINKS/SAMPLE_LINK/XREF_LINK/ID')
        for info in sample_info[:4]:
            #print info.text
            local.append(info.text)
        
        sample_attrib = root.findall('./SAMPLE/SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE/VALUE')
        for attrib in sample_attrib[:7]:
            #print attrib.text
            local.append(attrib.text)
        
        '''
        for country in root.findall('SAMPLE'):
            name = country.get('accession')
            print name
        '''
        print local
        resul.append('\t'.join(local))
    return resul


if to_do == "get":
    patterns = read_file(sys.argv[1])
    get_files(patterns)
    print "Done downloading xml"
elif to_do == "parse":
    files = glob.glob(sys.argv[1])
    reslist = parse_files(files)
    print "Done parsing xml"
    writeLIST_to_file(reslist, "metadata.txt")