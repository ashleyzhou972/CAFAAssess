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
# Change: use goatools (Haibo Tang) instead of Ontology
# Because goatools addressess alt_id
# Need to use my fork, for now, since goatools does not do edge-types, and the original goatools
# only uses is_a relationships
#
#
import os, sys
import argparse
import urllib
from collections import defaultdict
#from Ontology.IO.OboIO import OboReader
from goatools.obo_parser import GODag
import tempfile


def read_ontology(obo_path,verbose=False):
    if verbose:
        sys.stderr.write("Reading OBO file %s..." % obo_path)
    ontology = GODag(obo_path)
    return ontology

def get_ancestor_ids(ontology,alt_ids=False,verbose=False):
    # Returns ancestors IDs and alt_ids.
    # So this should not be used to actually count ancestors
    ancestors = defaultdict(set)
    for term_id in ontology:
        for parent in ontology[term_id].parents:
            ancestors[term_id].add(parent.id)
            if alt_ids: # Include alternative ids 
                ancestors[term_id] |= set(parent.alt_ids)
    return ancestors
 
def write_ancestors(ancestors, outpath):
    with open(outpath,"w") as outfile:
        for term_id in ancestors:
            outfile.write( "%s\t%s\n" % (term_id, ",".join(ancestors[term_id])) )

# Todo: write an ancestors + edges file using something like this loop
#>>> for i in list(ont.get_ancestors('GO:0009657')):
#...     for j in list(ont.get_node(i).succ):
#...         print(i,j.data,j.to_node.label)
#... 
# 
def write_ancestors_with_edges(ancestors, outpath):
    pass

def main(args):
    ontology = read_ontology(args.obopath, args.verbose)
    ancestors = get_ancestor_ids(ontology, args.verbose)
    write_ancestors(ancestors, args.outpath)
    

if __name__ == '__main__':
    # create the argument parser for this CLI
    parser = argparse.ArgumentParser(prog='obo2ancestors2',
        description='Parses an obo file and creates a file listing the ancestors of each accession ID')

    # setup the arguments for the argparser
    parser.add_argument('-i', '--input', required=True, dest='obopath', help='The obo file to parse')
    parser.add_argument('-o', '--output', required=True, dest='outpath', help='Output file')
    parser.add_argument('-v', '--verbose', action='store_true', default=False)
    parser.add_argument('-e', '--edges', default=False, help='Add edge type to each ancestor')

    # parse the arguments for this CLI
    args = parser.parse_args()
    main(args)

