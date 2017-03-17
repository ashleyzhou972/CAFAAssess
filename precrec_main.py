# -*- coding: utf-8 -*-
"""
Command-line interface of the precision-recall calculations in
CAFA assessment tool
@author: Ashley Zhou
last updated 03/17/2017
"""

import argparse
import tests
from precrec.precRec import PrecREC,read_benchmark
from precrec.GOPred import GOPred
import matplotlib.pyplot as plt
import numpy
import os



def get_namespace_index(namespace):
    '''
    convert namespace into indices
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
        print(namespace)
    return num


def taxon_name_converter(taxonID):
    #convert from taxonomy ID to name (i.e. from 9606 to HUMAN）
    taxonTable = {'10116':'RAT','9606':'HUMAN','3702':'ARATH','7955':'DANRE','44689':'DICDI','7227':'DROME','83333':'ECOLI','10090':'MOUSE','208963':'PSEAE','4896':'SCHPO','4932':'YEAST'}
    return taxonTable[taxonID]    
  
if __name__=='__main__':
    
    parser = argparse.ArgumentParser(description='Precision- Recall assessment for CAFA predictions.', )
    parser.add_argument('file',type=open,
                        help='Input the path of the prediction file. Filename should follow CAFA formats')
    #CAFA3 raw submission filename formats are listed here:https://www.synapse.org/#!Synapse:syn5840147/wiki/402192
    #example filename format: Doegroup_1_9606.txt/Doegroup_2_hpo.txt
    #If prediction file is already split by ontology it should follow Doegroup_1_9606_BPO.txt(or _MFO, _CCO)                  
    
    #parser.add_argument('plotfile', help='Input path+filename to save the PR plot')

    parser.add_argument('type',help = 'Input evaluation type: No Knowledge or Limited Knowledge', choices=['type1','type2','all'])    
    parser.add_argument('teamnum',help='Input team number',type=int)
    parser.add_argument('mode', help = 'Input the evaluation mode: full or partial', choices = ['full','partial'])
    args = parser.parse_args()
    
    print('Evaluating %s.\n' % args.file.name)
    
    #first split the prediction file into three ontologies
    all_pred = GOPred()
    pred_path = args.file
    obo_path = './precrec/go_20130615-termdb.obo'
    all_pred.read_and_split_and_write(obo_path,pred_path)
    info = [all_pred.author,all_pred.model,all_pred.keywords,all_pred.taxon]
    print('author: %s\n' % info[0])
    print('model: %s\n' % info[1])
    print('keywords: %s\n' % info[2][0])
    print('species:%s\n' % info[3])
    
    #parse file name 
    namefields = os.path.basename(pred_path.name).split('.')[0].split('_')
        
    #read benchmark for three ontologies
    for onto in ['bpo','cco','mfo']:
        print('ontology: %s\n' % onto)
        b = read_benchmark(onto, taxon_name_converter(namefields[2]),args.type,'./precrec/benchmark/',obo_path)
        path = os.path.splitext(pred_path.name)[0]+'_'+onto.upper()+'.txt'
        c = PrecREC(b,path)
        if c.exist:
            fm = c.Fmax_output(args.mode)
            print('fmax: %s\n' % fm[2])
            print('threshold giving fmax: %s\n' % fm[3])
            print('coverage: %s\n' % fm[4])
            yx = tests.read_Cafa2_sheet(onto, info[3],args.teamnum,info[1],args.type,args.mode)
            if yx[1]==round(fm[2],3):
                print('Fmax Match!\n')
                if yx[0]==round(fm[4],2):
                    print('Coverage Match!\n')
                    if yx[2]==fm[3]:
                        print('Threshold Match!\n')
