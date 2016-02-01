#!/usr/bin/env python
#
# Command line interface for the obo2ancestors script that allows you to take an
# obo file and list the ancestors of every term within the file.
#
# uses the Ontology module from Kamil Koziara
# @author Iddo Friedberg
# Tests run:
# 1. Separately tested using gawk script that it reads all terms in go-basic.obo
# 2. Spot-checked that all ancestors are recalled.
#
import os, sys
import argparse
import urllib
from collections import defaultdict
from Ontology.IO.OboIO import OboReader
import tempfile


def read_ontology(obo_path,verbose=False):
    ancestors = defaultdict(list)
    if verbose:
        sys.stderr.write("Reading OBO file %s..." % obo_path)
    obo_parser = OboReader(open(obo_path))
    ontology = obo_parser.read()
    return ontology

def get_ancestors(ontology,verbose=False):
    ancestors = {}
    if verbose:
        sys.stderr.write("Extracting ancestors....")
    for node in ontology.nodes:
        ancestors[node] = ontology.get_ancestors(node)
    return ancestors

# todo
def get_ancestors_with_edges(ontology, verbose=False):
    for node in ontology.nodes:
        for i in list(ontology.get_ancestors(node.label)):
            for j in list(ontology.get_node(i)):
                pass #TODO
def write_ancestors(ancestors, outpath):
    with open(outpath,"w") as outfile:
        for go, ancestors in list( ancestors.items() ):
            outfile.write("%s\t%s\n" % (go, ",".join(ancestors)) )

# Todo: write an ancestors + edges file using something like this loop
#>>> for i in list(ont.get_ancestors('GO:0009657')):
#...     for j in list(ont.get_node(i).succ):
#...         print(i,j.data,j.to_node.label)
#... 
# 
def write_ancestors_with_edges(ancestors, outpath):
    with open(outpath,"w") as outfile:
        for go, ancestors in list( ancestors.items() ):
            
            outfile.write("%s\t%s\n" % (go, ",".join(ancestors)) )

def main(args):
    ontology = read_ontology(args.obopath, args.verbose)
    ancestors = get_ancestors(ontology, args.verbose)
    write_ancestors(ancestors, args.outpath)
    

if __name__ == '__main__':
    # create the argument parser for this CLI
    parser = argparse.ArgumentParser(prog='obo2ancestors',
        description='Parses an obo file and creates a file listing the ancestors of each accession ID')

    # setup the arguments for the argparser
    parser.add_argument('-i', '--input', required=True, dest='obopath', help='The obo file to parse')
    parser.add_argument('-o', '--output', required=True, dest='outpath', help='Output file')
    parser.add_argument('-v', '--verbose', action='store_true', default=False)
    parser.add_argument('-e', '--edges', default=False, help='Add edge type to each ancestor')

    # parse the arguments for this CLI
    args = parser.parse_args()
    main(args)

