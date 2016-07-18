# -*- coding: utf-8 -*-
"""
Created on Thu Jul 14 15:02:03 2016
This creates the command-line interface of the precision-recall calculations in
CAFA assessment tool
@author: Ashley Zhou
"""

import sys
import os
sys.path.append(os.getcwd())
import argparse
from CAFAAssess.precRec import PrecREC,read_benchmark
from CAFAAssess.precrec.GOPred import GOPred
import matplotlib.pyplot as plt
import numpy

'''
arguments supplied to main:
updated on 07/12/2016
    ontology: 'BPO', 'MFO', 'CCO', 'HPO'
    team number: e.g. 117
    species: use taxon ID e.g. 9606 for Homo Sapien
    model: which model of the submission, e.g. 1,2 or 3 
    saveplot path: the path where the PR curve plot is saved (should be an absolute path??)
    #use default go.obo: True or supply obopath (should be an absolute path??)
    #use default benchmark: True or supply benchmark path (Ontology-specific!!) (should be an absolute path)
    #submission folder: folder containing all submissions with each team in separate  folders titled by team number
'''

def get_namespace_index(namespace):
    '''
    copied from confidence.py 07/15/2016
    '''
    num = None
    if namespace=='BPO' or namespace=='bpo':
        num = 0
    elif namespace=='MFO' or namespace=='mfo':
        num = 1
    elif namespace=='CCO' or namespace=='cco':
        num =2
    else:
        raise ValueError("name space not found, check prediction files")
        print namespace
    return num
'''
def team_name_converter(team_number):
    
    #convert team number to team name for plotting purposes
    #07/18/2016: waiting for input from actual CAFA3
    
    return "teamName"
    
def taxon_name_converter(taxonID):
    return "taxon"
'''    

parser = argparse.ArgumentParser(description='Precision- Recall assessment for CAFA predictions.', )
parser.add_argument('ontology',help='Input ontology',choices=['BPO','MFO','CCO'])
parser.add_argument('team',help = 'Input team number',type=int)
parser.add_argument('taxon', help= 'Input taxon ID, this will only be used to name the plot', type=int)
#If it's all species combined, enter 0 as taxon ID(07/18/2016)
#parser.add_argument('model',help = 'Input model number', choices=['1','2','3'])
parser.add_argument('file',type=open,
                    help='Input the path of the prediction file. File should be split according to ontology, and should be a .txt file')
parser.add_argument('plotfile', help='Input path+filename to save the PR plot')
#parser.parse_args(['BPO','117','./CAFAAssess/confidence/117/PaccanaroLab_1_9606_BPO.txt','./CAFAAssess/confidence/117/'])
args = parser.parse_args()

print('Evaluating %s.\n' % args.file.name)
print('Ontology: %s\n' % args.ontology)
print('Species: %s\n' % args.taxon)
print('Team: %s\n' % args.team)
#print('model: %s\n' % args.taxon)

bench = [read_benchmark('BPO'), read_benchmark('MFO'), read_benchmark('CCO')]
index = get_namespace_index(args.ontology)
all_pred = GOPred()
all_pred.read(args.file)
c = PrecREC(bench[index],all_pred)
fm = c.Fmax_output(99)

plt.plot(fm[1],fm[0])
plt.axis([0,1,0,1])
plt.yticks(numpy.arange(0,1,0.1))
plt.xlabel('Recall')
plt.ylabel('Precision')
plt.title(str(args.team)+" "+str(args.taxon)+" "+args.ontology)
plt.savefig(args.plotfile,dpi=200)
plt.close()

print('fmax value for this prediction is: %s.\n' % fm[2])
print('PR plot is saved to %s.\n' % args.plotfile)
