# -*- coding: utf-8 -*-
"""
Command-line interface of the precision-recall calculations in
CAFA assessment tool
@author: Ashley Zhou
last updated 03/17/2017
TODO: positional argument numbers and help messages
type='all'
plot more than one curve
more than one prediction file input
"""

import argparse
from precrec.precRec import PrecREC,read_benchmark
from precrec.GOPred import GOPred
import numpy
import os
import matplotlib.pyplot as plt


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


def plotSingle(precision, recall, author, model, taxon, ontology, mode, TYPE):
    '''
    precision should be a list of precision values 
    recall should be a list of recall values
    pay attention to order of precision-recall
    recall is on x-axis, but supplied second
    '''
    plt.plot(recall,precision)
    plt.axis([0,1,0,1])
    plt.yticks(numpy.arange(0,1,0.1))
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    figuretitle = author+" "+str(model)+" "+str(taxon)+ " " + ontology+' '+ 'mode:'+mode+' '+'type:'+TYPE
    plt.title(figuretitle)
    figurename = './plots/'+author+"_"+str(model)+"_"+str(taxon)+ "_" + ontology+'_'+ mode+'_'+TYPE+'.png'
    plt.savefig(figurename,dpi=200)
    plt.close()

def plotMultiple(title,*args):
    for i in args:
        plt.plot(i[1],i[0])
    plt.axis([0,1,0,1])
    plt.yticks(numpy.arange(0,1,0.1))
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.legend([i[2] for i in args], loc='best')
    plt.title(title)
    if title == None:
        figurename = './plots/Combined_plot.png'
    else:
        figurename = './plots/'+title+'.png'
    
    plt.savefig(figurename,dpi=200)
    plt.close()
def typeConverter(oldType):
    if oldType=='type1':
        newType = 'No_Knowledge'
    elif oldType == 'type2':
        newType = 'Limited_Knowledge'
    elif oldType == 'all':
        newType = 'All'
    return(newType)


if __name__=='__main__':
    
    parser = argparse.ArgumentParser(description='Precision- Recall assessment for CAFA predictions.', )
    parser.add_argument('file',type=open,
                        help='Input prediction file. Filename should follow CAFA formats. Accepts more than one predictions.',
                        nargs = '+')
    #CAFA3 raw submission filename formats are listed here:https://www.synapse.org/#!Synapse:syn5840147/wiki/402192
    #example filename format: Doegroup_1_9606.txt/Doegroup_2_hpo.txt
    #If prediction file is already split by ontology it should follow Doegroup_1_9606_BPO.txt(or _MFO, _CCO)                  
    
    
    parser.add_argument('-t','--t',dest='type',help = 'Input evaluation type: No Knowledge or Limited Knowledge', choices=['type1','type2','all'],required=True)    
    parser.add_argument('-o','--o', dest= 'obo_path',help = 'Input the obo file path',default = './precrec/go_20130615-termdb.obo')
    parser.add_argument('-m','--m',dest='mode', help = 'Input the evaluation mode: full or partial', choices = ['full','partial'],required = True)
    parser.add_argument('-b','--b',dest='bfolder', help = 'Input the path to the benchmark folder, default CAFA2 benchmarks provided', default = './precrec/benchmark/')
    parser.add_argument('-title', dest='title',help = 'Input title of combined plot, if multiple prediction files are supplied',default = ' ')
    args = parser.parse_args()
    os.mkdir('./plots/')
    os.mkdir('./results/')
    
    num = len(args.file)
    for f in args.file:
        print('Evaluating %s.\n' % f.name)
        result = open("./results/%s_results.txt" % os.path.basename(f.name)[0],'w')
        #first split the prediction file into three ontologies
        all_pred = GOPred()
        pred_path = f
        obo_path = args.obo_path
        benchmarkFolder = args.bfolder
        all_pred.read_and_split_and_write(obo_path,pred_path)
        info = [all_pred.author,all_pred.model,all_pred.keywords,all_pred.taxon]
        print('AUTHOR: %s\n' % info[0])
        result.write('AUTHOR:%s\n' % info[0])
        print('MODEL: %s\n' % info[1])
        result.write('MODEL: %s\n' % info[1])
        print('KEYWORDS: %s\n' % info[2][0])
        result.write('KEYWORDS: %s\n' % info[2][0])
        print('Species:%s\n' % info[3])
        result.write('Species:%s\n' % info[3])
        print('benchmark type:%s\n' % typeConverter(args.type)) 
        result.write('benchmark type:%s\n' % typeConverter(args.type))
        print('mode:%s\n' % args.mode)
        result.write('mode:%s\n' % args.mode)
        result.write('%s:\t%s\t%s\t%s\n' % ('Ontology','Fmax','Threshold','Coverage'))
        #parse file name 
        namefields = os.path.basename(pred_path.name).split('.')[0].split('_')
            
        #read benchmark for three ontologies
        for onto in ['bpo','cco','mfo']:
            print('ontology: %s\n' % onto)
            b = read_benchmark(onto, taxon_name_converter(namefields[2]),args.type,benchmarkFolder,obo_path)
            path = os.path.splitext(pred_path.name)[0]+'_'+onto.upper()+'.txt'
            c = PrecREC(b,path)
            if c.exist:
                fm = c.Fmax_output(args.mode)
                print('fmax: %s\n' % fm[2])
                print('threshold giving fmax: %s\n' % fm[3])
                print('coverage: %s\n' % fm[4])
                plotSingle(fm[0],fm[1],info[0],info[1],info[3],onto,args.mode,typeConverter(args.type))
                result.write('%s:\t%s\t%s\t%s\n' % (onto,fm[2],fm[3],fm[4]))
        result.close()
