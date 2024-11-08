#!/usr/bin/env python2.7
# -*- coding: UTF-8 -*-
# Function:
# 过滤run_gatk_pileup_for_sample.py结果，作为verify_concordance.py的输入
# Version: 0.1
# Email: zhengshimao007@163.com
# Author: Shimao Zheng


import optparse
import sys
import os

desc = """Program to filter GATK Pileup on a single sample"""
parser = optparse.OptionParser(version='%prog version 0.1 08/November/2024', description=desc)
parser.add_option('-i', '--pileup', help='pileup file [ required ]', action='store')
parser.add_option('-o', '--outfile', help='filtered pileup file for script verify_concordance.py [ <pileup>4concordance ]', action='store')
parser.add_option('-e', '--emptyfile', help='empty pileup file [ <pileup>.empty ]', action='store')

(opts, args) = parser.parse_args()

pileup = opts.pileup
if not pileup :
    parser.print_help()
    sys.exit(1)

if not os.path.exists(pileup):
    print('ERROR: Input pileup file {0} cannot be found.'.format(pileup))
    sys.exit(1)


outfile = opts.outfile
if not outfile :
	outfile = pileup + "4concordance"

emptyfile = opts.emptyfile
if not emptyfile:
	emptyfile = pileup + ".empty"

with open(pileup, 'r') as f1, \
	 open(outfile,'w') as f2, \
	 open(emptyfile,'w') as f3:
	for line in f1:
		columns = line.split(' ')
		if columns[4] != "":
			f2.write(line) # not empty
		else:
			f3.write(line) # empty
