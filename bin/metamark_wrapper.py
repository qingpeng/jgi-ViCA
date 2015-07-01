#!/usr/bin/env python

import sys
import simplejson as json
import argparse
import os
import subprocess

parser = argparse.ArgumentParser(description='A script to generate a feature matrix  \
using emmission data from Metamark')
parser.add_argument('-r','--root', help='a root input directory', required=True)
parser.add_argument('-p','--par', help='an alternate location for the metamark parameter files', default=os.path.join(os.path.dirname(os.path.abspath(sys.argv[0])),"gm_parameters"))
parser.add_argument('-c,','--config',help= 'a json formated config file', default='config.json')
args = parser.parse_args()
