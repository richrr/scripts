import os
import sys
from utils import *
import re

infile = sys.argv[1]
lines = read_file(infile)

def escapeDoubleUnderscores(lines):
    esc_list = list()
    for l in lines:
        l = l.strip()
        esc_l = l.replace('__','\_\_')
        esc_list.append(esc_l)
    return esc_list

def checkWeirdColon(l):
    if ':' in l and re.search(':\S+?:',l) : # search any colon NOT followed by white space
        return True
    return False

def fix_WeirdColon(l):
    l = re.sub(r':(\d)', r'-\1', l)  # r'' is always a safe choice. r'' prevents converting \s and \n to space and newline, but \d is fine
    l = re.sub(r'(\\_\\_):(\w+):', r'\1\\\\[\2]', l)  #\_\_:Chthoniobacterales: becomes \_\_\\[Chthoniobacterales]
    return l

escaped_lines = escapeDoubleUnderscores(lines)

taxa_to_check = list()
fixed_taxa = list()

for l in escaped_lines:
    if "Other" in l:
        taxa_to_check.append(l)
        continue
    if checkWeirdColon(l):
        taxa_to_check.append(l)
        continue
    if "__" in l:
        taxa_to_check.append(l)
        continue
    if ';' in l:
        taxa_to_check.append(l)
        continue
    fixed_taxa.append(l)

print "Check the following taxa" , '\n\t', '\n\t'.join(taxa_to_check)

for l in taxa_to_check:
    l = fix_WeirdColon(l)
    l = l.replace(":Other",'')
    l = l.replace(" :",': ')
    if 'g\_\_' not in l:
        l = l + '$'
    fixed_taxa.append(l)

print "Use the following fixed taxa instead" , '\n\t', '\n\t'.join(fixed_taxa)

writeLIST_to_file(fixed_taxa, infile+'-fixed-taxa.txt' )



