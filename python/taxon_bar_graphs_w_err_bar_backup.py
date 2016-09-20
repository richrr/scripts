import os
import sys
sys.path.insert(0,"/usr/lib/python2.7/dist-packages/")
from utils import *
import operator
from time import localtime, strftime
import argparse
import math
import numpy
import matplotlib.pyplot as plt
import brewer2mpl
import pylab
from scipy import stats
#from statsmodels.stats import multicomp

#usage: python ~/scripts/python/taxon_bar_graphs_w_err_bar.py -i table_mc80000_sorted_L2.txt -g ../../../../combined_rosana_qiime_mapping_file.txt -c Group
# python ~/scripts/python/taxon_bar_graphs_w_err_bar.py -i sorted_row_sum_60samples-013113MWITS-pr.fa.otus.faOTU_even250_L2.txt -g ../../../../../its-michigan-mapping-file.txt -c Season SoilAgeYears
#python ~/scripts/python/taxon_bar_graphs_w_err_bar.py -i otu_table_mc2_w_tax_even3400_L3.txt -g ../../../sample_sheet.txt -c Group -o L3_ -v unequal
#python ~/scripts/python/taxon_bar_graphs_w_err_bar.py -i otu_table_mc2_w_tax_even3400_L3.txt -g ../../../sample_sheet.txt -c Group -o L3_ -s manwu


def main(args):

    parser = argparse.ArgumentParser(description='Taxon summary plots with error bars')
    parser.add_argument('-i', '--infile')
    parser.add_argument('-o', '--outfile', default="outfig")  # output filename
    parser.add_argument('-c', '--categories' , nargs='*') \
        # creates a list of groups (categories) to be summarized (enter as space delimited on cmd line), currently only one category
    parser.add_argument('-g', '--groupsfile')  # file containing the group where each sample belongs
    parser.add_argument('-r', '--contrastsfile', default='contrasts.txt')  # 2 col. tab-delim file containing the 2 groups to be compared
    parser.add_argument('-f', '--pvaladjmethod', default='bf')  # pval adjust method: bonferroni (bf), benjamini hochberg (bh) 
    parser.add_argument('-v', '--variance', default='equal')  # equal or unequal variance
    parser.add_argument('-d', '--delimiter', default='\t')  # delimiter for file
    parser.add_argument('-t', '--imagetype', default='pdf')  # file type of images, default pdf, other allowed are png or jpg
    parser.add_argument('-s', '--statstest', default='ttest')  # mann whitney u test (manwu), or t test (ttest)

    args = parser.parse_args()

    if len(sys.argv)==1 :
        parser.print_help()
        sys.exit('\natleast one argument required\n')

    infile  = args.infile
    outpfile = args.outfile
    delim = args.delimiter
    f_type = args.imagetype
    contrasts = read_file(args.contrastsfile)
    variance = args.variance
    p_val_adj = args.pvaladjmethod
    stats_test = args.statstest
    categories_to_print = ['Group']
    if args.categories != None:
        categories_to_print = args.categories
    
    global logfile, logfilelist
    logfile = 'tax-plot-err-bar-log-file.txt'
    logfilelist = list()

    group_file = ''
    if args.groupsfile != None:
        group_file = args.groupsfile
    mapping_file = read_file(group_file) 
    
    lines = read_file(infile)
    if 'Taxon' not in lines[0]:
        sys.exit('\nheader line is required, tab-delimited, starts with taxon, followed by samples\n')

    for c in categories_to_print:
        draw_tax_plots_w_err(c, mapping_file, lines, infile, outpfile, delim, f_type, contrasts, p_val_adj, variance, stats_test)

def draw_tax_plots_w_err(c, mapping_file, lines, infile, outpfile, delim, f_type, contrasts, p_val_adj, variance, stats_test):
    # find index of SampleID and category to be summarized
    labels = mapping_file[0].split(delim)
    indx_sampleid = indx_categ = ''

    try:
        indx_sampleid = labels.index("#SampleID")
    except ValueError:
        print "SampleID not in list."

    try:
        indx_categ = labels.index(c)
        #print indx_categ
    except ValueError:
        print "%s not in list." %c

    categ_samples_dict = list_to_dict(mapping_file, delim, ',', "current", indx_categ, indx_sampleid)
    #print categ_samples_dict
    category_list_sorted = sorted(categ_samples_dict) # so the categories are sorted (alphabetically)
    #category_list_sorted = ["0", "Early", "Middle", "Late"]
    #category_list_sorted = ["0", "105", "155", "210", "450", "845", "1410", "2510", "3210", "4010"]
    #print categ_samples_dict, categ_samples_dict_sorted
    
    # collect samples belonging to the same group
    categ_col_indx_dict = dict() # 0 based indexes, key is category, value is cols of samples in that category
    taxon_values_mean = dict() # 0 based indexes, key is taxon, value is tuple of means of values of samples in same category i.e. (mean_categ1, mean_categ2, mean_categ3)
    taxon_values_stderr = dict() # 0 based indexes, key is taxon, value is tuple of std err of values of samples in same category i.e. (err_categ1, err_categ2, err_categ3)
    for l in lines:
        #print l, "------------------"
        l = l.strip()
        contents = l.split(delim)
        if '#' in l or not l:
            continue
        if 'Taxon' in l:
            for k in category_list_sorted:
               v = categ_samples_dict[k]
               indxs = [contents.index(i) for i in v.split(',') if i in contents]
               categ_col_indx_dict[k] = indxs
            continue
        taxon_name = contents[0]
        taxon_values_mean[taxon_name] = list()
        taxon_values_stderr[taxon_name] = list()
        group_values_dict = dict()
        # calc. mean, std. err for each group
        for k in category_list_sorted:
            v = categ_col_indx_dict[k]
            vals = [float(contents[i]) for i in v]
            group_values_dict[k] = vals
            #print taxon_name, k, v, vals
            taxon_values_mean[taxon_name].append(str(numpy.mean(vals)))  # convert these directly to float
            taxon_values_stderr[taxon_name].append(str(numpy.std(vals)/math.sqrt(len(vals))))
        #print taxon_values_mean, taxon_values_stderr
        two_sample_stats_test(taxon_name, group_values_dict, contrasts, p_val_adj, variance, delim, len(lines)-1, stats_test) # the first line is header
    #print taxon_values_mean, taxon_values_stderr
    draw_stacked_bar_w_err_bar(c, taxon_values_mean, taxon_values_stderr, category_list_sorted, outpfile, f_type)


def two_sample_stats_test(taxon_name, group_values_dict, contrasts, p_val_adj, variance, delim, numb_taxa, stats_test):
    for c in contrasts:
        c = c.strip()
        if '#' in c or not c:
            continue
        categ1, categ2 = c.split(delim) # categories to be compared
        #print categ1, categ2
        if categ1 not in group_values_dict or categ2 not in group_values_dict:
            continue
        group1 = group_values_dict[categ1]
        group2 = group_values_dict[categ2]
        #print group1, group2
        if stats_test == 'manwu':
            two_sample_manwu = ''
            try:
                two_sample_manwu = stats.mannwhitneyu(group1, group2)
            except: #ValueError: All numbers are identical in amannwhitneyu
                continue
            mstats, pval_ = two_sample_manwu
            if p_val_adj == 'bf':
                pval_adjusted_c = pval_ * len(contrasts)
                pval_adjusted_t = pval_ * numb_taxa
                pval_adjusted = pval_ * len(contrasts) * numb_taxa
                if pval_ <= 0.05:#unadjusted
                    print "%s\t%s\tmanwu-stat\t%.3f\t p-val\t%.3f\tbf adj p-val(contrasts)\t%.3f\tbf adj p-val(taxa)\t%.3f\tbf adj p-val(contrasts & taxa)\t%.3f" % (c, taxon_name, mstats, pval_, pval_adjusted_c, pval_adjusted_t, pval_adjusted)
            continue
        #return value of t test is a tuple containing the t-statistic and the p-value and these are the results of a two-sided test.
        if variance == 'unequal':
            # assuming unequal population variances
            two_sample_diff_var = stats.ttest_ind(group1, group2, equal_var=False)
            tstats, pval_ = two_sample_diff_var
            if p_val_adj == 'bf':
                pval_adjusted_c = pval_ * len(contrasts)
                pval_adjusted_t = pval_ * numb_taxa
                pval_adjusted = pval_ * len(contrasts) * numb_taxa
                if pval_ <= 0.05:#unadjusted
                    print "%s\t%s\tt-stat\t%.3f\t p-val\t%.3f\tbf adj p-val(contrasts)\t%.3f\tbf adj p-val(taxa)\t%.3f\tbf adj p-val(contrasts & taxa)\t%.3f" % (c, taxon_name, tstats, pval_, pval_adjusted_c, pval_adjusted_t, pval_adjusted)
            continue
        two_sample = stats.ttest_ind(group1, group2, equal_var=True)
        tstats, pval_ = two_sample
        if p_val_adj == 'bf':
            pval_adjusted_c = pval_ * len(contrasts)
            pval_adjusted_t = pval_ * numb_taxa
            pval_adjusted = pval_ * len(contrasts) * numb_taxa
            if pval_ <= 0.05:
                print "%s\t%s\tt-stat\t%.3f\t p-val\t%.3f\tbf adj p-val(contrasts)\t%.3f\tbf adj p-val(taxa)\t%.3f\tbf adj p-val(contrasts & taxa)\t%.3f" % (c, taxon_name, tstats, pval_, pval_adjusted_c, pval_adjusted_t, pval_adjusted)

#------ TO DO ----------------
# do for multiple input files
# sort and put the most dominant on top
# most dominant 5-10 only
def draw_stacked_bar_w_err_bar(c, means_d, stderr_d, category_labels, outpfile, f_type):
    N = len(category_labels)
    # Get colors from ColorBrewer (all colorbrewer scales: http://bl.ocks.org/mbostock/5577023)
    custom_color_list = list()
    custom_color_list.extend(brewer2mpl.get_map('Paired', 'qualitative', 10).mpl_colors)
    custom_color_list.extend(brewer2mpl.get_map('Set2', 'qualitative', 6).mpl_colors)
    custom_color_list.extend(brewer2mpl.get_map('Set1', 'qualitative', 7).mpl_colors)
    custom_color_list.extend(brewer2mpl.get_map('Set3', 'qualitative', 10).mpl_colors)
    
    # Set the random seed for consistency
    numpy.random.seed(12)
    ind = numpy.arange(N)    # the x locations for the groups
    width = 0.35       # the width of the bars
    
    bottom_value = [0 for m in category_labels]
    legend_bar = list()
    legend_text = list()
    counter = 0
    #taxons_sorted_list = sorted(means_d) # this sorts the taxons (keys of the means_d dictionary) and saves it as a list

    # print in bar graph (bottom to top), the taxons in decreasing order of average abundance across groups.
    #http://www.quora.com/How-do-I-sum-the-values-of-each-key-in-a-Python-dictionary-and-then-display-the-result-as-key-summed-value
    d_total_list_values = dict((k, sum([float(i) for i in v])/float(len(v))) for k,v in means_d.items()) # convert values in list to float, avg them, assign total as value to key
    taxons_sorted_list = sorted(d_total_list_values, key=d_total_list_values.get, reverse=True)# this sorts the taxons (keys) as per values and saves it as a list
    #print d_total_list_values, taxons_sorted_list
    for t in taxons_sorted_list:
        if counter == len(custom_color_list):  # cycle through the list of available colors if we have more taxons
            counter = 0
        mean = [float(i) for i in means_d[t]]
        stderr = [float(i) for i in stderr_d[t]]
        p = plt.bar(ind, mean, width, color=custom_color_list[counter], yerr=stderr, bottom=bottom_value)
        bottom_value = [x+y for x,y in zip(bottom_value, mean)]
        legend_bar.append(p[0])
        legend_text.append(t)
        counter += 1
    plt.ylabel('Cumulative relative abundance')
    plt.xlabel(c)
    plt.title('Taxonomic summary')
    plt.xticks(ind+width/2., tuple(category_labels) )
    #plt.legend( tuple(legend_bar), tuple(legend_text) )
    figname = outpfile + '_' + c + '.' + f_type
    plt.savefig(figname, close=True)

    # Klaus se at stackexchange
    #http://stackoverflow.com/questions/4534480/get-legend-as-a-seperate-picture-in-matplotlib
    legend_fig = pylab.figure()
    legend = pylab.figlegend( tuple(legend_bar), tuple(legend_text), loc = 'center' )
    legend_fig.canvas.draw()
    legend_figname = 'legend_' + figname 
    legend_fig.savefig(legend_figname, bbox_inches=legend.get_window_extent().transformed(legend_fig.dpi_scale_trans.inverted()))
    pylab.close()



if __name__=='__main__':
    datetime = strftime("%a, %d %b %Y %I:%M:%S %p", localtime())
    cmd = 'echo ' + datetime 
    os.system(cmd)
    main(sys.argv)

