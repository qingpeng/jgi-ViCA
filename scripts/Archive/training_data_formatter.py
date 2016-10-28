#!/usr/bin/env python

import argparse
import simplejson as json
import os
import subprocess
import numpy


parser = argparse.ArgumentParser( description='A script to select taxonomic groups with particular attributes from a reftree feature vector database and split the data into test and train datasets' )
parser.add_argument('-c', '--config', help ="A json formatted config file")
parser.add_argument('-t', '--trainfile', help ="An output vector file", type=argparse.FileType('w'), default='-')
parser.add_argument('-e', '--testfile', help ="An output vector file", type=argparse.FileType('w'), default='testdata.txt')
parser.add_argument('-s', '--switch', help ="reftree or vector to use") 
parser.add_argument('-g', '--segment', help ="segment or genome, based vector") 
parser.add_argument('-r', '--reftree', help ="Reftree database directory location") 
parser.add_argument('-v', '--vector', help ="Vector file if no reftree database") 
args = parser.parse_args()

# Read the configuration file
config = json.load(open(args.config, 'r'))["training_data_formatter"]

# Functions

def validate(a,b):
    """A function to validate that every item in dictionary 'a' is present in 'b' and 
    its value matched the value in 'a' including cases where that value is 'None in 'a'
     but not present in 'b'. Returns True or False."""
    for key, val in a.iteritems():
        if val == "None":
            if key in b and b[key] != "None":
                return False
        else:
            if key not in b or b[key] != val:
                return False
    return True 
        
def linetodict_vector(line): # based on genome
    print line
    """split the data vector into a dictionary for vector file directly"""
    lv =line.strip().split("\t")
    if len(lv) < 5:
        return None
    lv1_fields = lv[1].split(" ")
    print lv1_fields
    if len(lv1_fields) <2:
        return None
    attributes = lv1_fields[1]
#    c3 = lv[2].split(" ")
#    attributes = lv[1]+","+ c3[2]
    attlist = attributes.split(',')
    print attlist
    d = dict(s.split('=') for s in attlist)
    taxlist = (d["taxonomy"].split('/'))
    t = dict(r.split(':') for r in taxlist)
    d["taxid"] = lv[0]
    d["taxonomy"] = t
#    d["readId"] = c3[0]
    d["refseqId"] = lv1_fields[0]
    d["vector"] = lv[2:]
    return d
    
def linetodict_vector_segment(line): #based on segments
#    print line
    """split the data vector into a dictionary for vector file directly"""
    lv =line.strip().split("\t")
    if len(lv) < 5:
        return None
    lv1_fields = lv[1].split(" ")
#    print lv1_fields
    if len(lv1_fields) <3:
        return None
    attributes = lv1_fields[2]
#    c3 = lv[2].split(" ")
#    attributes = lv[1]+","+ c3[2]
    attlist = attributes.split(',')
#    print attlist
    d = dict(s.split('=') for s in attlist)
    taxlist = (d["taxonomy"].split('/'))
    t = dict(r.split(':') for r in taxlist)
    d["taxid"] = lv[0]
    d["taxonomy"] = t
#    d["readId"] = c3[0]
    d["refseqId"] = lv1_fields[0]
    d["vector"] = lv[2:]
    return d
    

def linetodict(line):
    """split the data vector into a dictionary"""
    lv =line.strip().split("\t")
    if len(lv) < 5:
        return None
    c3 = lv[2].split(" ")
    attributes = lv[1]+","+ c3[2]
    attlist = attributes.split(',')
    d = dict(s.split('=') for s in attlist)
    taxlist = (d["taxonomy"].split('/'))
    t = dict(r.split(':') for r in taxlist)
    d["taxid"] = lv[0]
    d["taxonomy"] = t
    d["readId"] = c3[0]
    d["refseqId"] = c3[1]
    d["vector"] = lv[3:]
    return d
    

def filter(d, params):
    """Returns a class label and vector given  dictionaries of params and data"""
    # create a dictionary of lists containing the unique parameters for each taxid (needed to define wanted_keys)
    taxcatdict = {}
    for key1, val1 in params.iteritems(): #  "ssRNAPhage" : {"taxcat" : "439488", "gencode" : "11"},
        tmplist =[]
        for key2 in val1:
            tmplist.append(key2) # "taxcat", "gencode"
        tid = params[key1]["taxcat"] # "439488"
        if tid in taxcatdict:
            taxcatdict[tid] = list(set(taxcatdict[tid] + tmplist))
        else:
            taxcatdict[tid] = tmplist 
    #iterate over the list pf parameters finding the category to assign the vector to   
    for key, val in params.iteritems():
        dc = d.copy()
        taxonomy = d["taxonomy"]
        if val["taxcat"] in taxonomy:
            dc["taxcat"] = val["taxcat"]
            if validate(val,dc) == True:
                dsub = {}
                wanted_keys = taxcatdict[val["taxcat"]] # The keys you want
                for item in wanted_keys:
                    if item in dc:
                        dsub[item] = dc[item]
                if validate(dsub, val) == True:
                    return [key] + dc["vector"]
                else:
                    continue
            else:
                continue
        else:
            continue
    

# MAIN

# validate config categories


if args.switch == "reftree":


    # Variables
    dbdir,dbname = os.path.split(os.path.abspath(args.reftree))

    # retrieve all data
    reftreeopts = ["reftree.pl", "--dbDir", dbdir , "--db", dbname , "--subtree" ,"1"] 
    p0 = subprocess.Popen(reftreeopts, stdout=subprocess.PIPE)

    input_stream = p0.stdout
    #

    # Create lists to hold taxids
    testlist = []
    trainlist = []

    # For each line 
    for line in input_stream:
        # Convert the line to a dictionary
        d = linetodict(line)
        # decide whether data should go in the test or train bin
        if d:
            if d["taxid"] in trainlist:
                out = args.trainfile
            elif d["taxid"] in testlist:
                out = args.testfile
            else:
                selection = numpy.random.choice(["train","test"], p=[1-float(config["testprop"]),float(config["testprop"])])
                if selection == "train":
                    out = args.trainfile
                    trainlist.append(d["taxid"])
                elif selection == "test":
                    out = args.testfile
                    testlist.append(d["taxid"])
                else:
                    raise ValueError('numpy.random.choice could not select a file to write the vector to')
            # Create a list with the category and vector if a match exists
            vectlist = filter(d, config["categories"])
            if vectlist: # If the vector list exists write it to the file
                out.write('\t'.join(vectlist))
                out.write('\n')
        else:
            continue
    
else:
    input_stream = open(args.vector,'r')

    # Create lists to hold taxids
    testlist = []
    trainlist = []

    # For each line 
    for line in input_stream:
        # Convert the line to a dictionary
        if args.segment == "genome":
            d = linetodict_vector(line)
        else:
            d = linetodict_vector_segment(line)
        # decide whether data should go in the test or train bin
        if d:
            if d["taxid"] in trainlist:
                out = args.trainfile
            elif d["taxid"] in testlist:
                out = args.testfile
            else:
                selection = numpy.random.choice(["train","test"], p=[1-float(config["testprop"]),float(config["testprop"])])
                if selection == "train":
                    out = args.trainfile
                    trainlist.append(d["taxid"])
                elif selection == "test":
                    out = args.testfile
                    testlist.append(d["taxid"])
                else:
                    raise ValueError('numpy.random.choice could not select a file to write the vector to')
            # Create a list with the category and vector if a match exists
            vectlist = filter(d, config["categories"])
            if vectlist: # If the vector list exists write it to the file
                out.write('\t'.join(vectlist))
                out.write('\n')
        else:
            continue
    