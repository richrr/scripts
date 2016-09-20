#!/usr/bin/env python
# File created on 09 Jul 2013
from __future__ import division
#https://groups.google.com/forum/#!searchin/qiime-forum/PCoA/qiime-forum/zigFP_wKaps/QPzKZgQW55IJ
# adapted from
__author__ = "Yoshiki Vazquez Baeza"
__copyright__ = "Copyright 2013, ApocaQIIME"
__credits__ = ["Yoshiki Vazquez Baeza"]
__license__ = "GPL"
__version__ = "1.7.0"
__maintainer__ = "Yoshiki Vazquez Baeza"
__email__ = "yoshiki89@gmail.com"
__status__ = "Use at your own risk"

# the creation of the convex hull is based on the example as provided by scipy
# see this URL http://docs.scipy.org/doc/scipy-dev/reference/generated/scipy.spatial.ConvexHull.html
import os
import sys
sys.path.insert(0,"/usr/lib/python2.7/dist-packages/")
from scipy.spatial import ConvexHull
import numpy as np
import matplotlib.pyplot as plt
from qiime.sort import natsort
from qiime.parse import parse_coords, parse_mapping_file
from qiime.util import (parse_command_line_parameters, make_option, qiime_system_call)
from qiime.colors import get_qiime_hex_string_color

script_info = {}
script_info['brief_description'] = ""
script_info['script_description'] = ""
script_info['script_usage'] = [("","","")]
script_info['output_description']= ""
script_info['required_options'] = [
    make_option('-i','--coordinates_fp', type="existing_filepath",help='the '
    'input coordinates filepath'),
    make_option('-m','--mapping_file_fp', type="existing_filepath",help='the '
    'mapping file filepath'),
    make_option('-c','--category',type="string",help='header name of the '
    'category of interest that you want the convex hulls to be created on')
]
script_info['optional_options'] = [
    make_option('-o','--output_fp', type="new_filepath", help='filename and '
    'format for the plot, the extension will determine the format of the '
    'output file (pdf, eps or png)', default='convex_hull.pdf')
]
script_info['version'] = __version__

def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    coordinates_fp = opts.coordinates_fp
    mapping_file_fp = opts.mapping_file_fp
    category_header_name = opts.category
    output_fp = opts.output_fp

    coords_headers, coords_data, coords_eigenvalues, coords_percents = parse_coords(open(coordinates_fp, 'U'))
    mapping_data, mapping_headers, _ = parse_mapping_file(open(mapping_file_fp, 'U'))

    category_header_index = mapping_headers.index(category_header_name)
    category_names = list(set([line[category_header_index] for line in mapping_data]))

    xtitle = 'PC1 (%.0f%%)' % round(coords_percents[0])
    ytitle = 'PC2 (%.0f%%)' % round(coords_percents[1])
    main_figure = plt.figure()
    main_axes = main_figure.add_subplot(1, 1, 1, axisbg='white')
    plt.xlabel(xtitle)
    plt.ylabel(ytitle)
    main_axes.tick_params(axis='y')
    main_axes.tick_params(axis='x')

    
    # sort the data!!! that way you can match make_3d_plots.py
    for index, category in enumerate(natsort(category_names)):
        sample_ids_list = [line[0] for line in mapping_data if line[category_header_index] == category]

        qiime_color = get_qiime_hex_string_color(index)

        if len(sample_ids_list) < 3:
            continue

        indices = [coords_headers.index(sample_id) for sample_id in sample_ids_list]
        points = coords_data[indices, :2]# * coords_percents[:2]

        hull = ConvexHull(points)
        main_axes.plot(points[:,0], points[:,1], 'o', color=qiime_color)
        for simplex in hull.simplices:
            main_axes.plot(points[simplex,0], points[simplex,1], 'k-')
        main_axes.plot(points[hull.vertices,0], points[hull.vertices,1], '--', lw=2, color=qiime_color)
        # plt.plot(points[hull.vertices[0],0], points[hull.vertices[0],1], '--', color=qiime_color)
    #plt.show()

    main_figure.savefig(output_fp)

if __name__ == "__main__":
    main()