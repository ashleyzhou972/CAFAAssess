# -*- coding: utf-8 -*-
"""
Created on Wed May 18 11:55:27 2016

@author: Ashley Zhou
"""


import sys
sys.path.append('/home/nzhou/git')
import os
from CAFAAssess.precRec import GOPred,benchmark,PrecREC
from Ontology.IO import OboIO
import gzip
import zipfile
import tarfile


pred_folder = "/home/nzhou/old computer/Documents/CAFA2/CAFA2_submissions/"  

#pred_path = open('/home/nzhou/git/CAFAAssess/precrec/M1HS.74.Homo_sapiens_BPO.txt')
#pred_path_ori = open('/home/nzhou/git/CAFAAssess/precrec/M1HS.74.Homo_sapiens.txt')



def read_benchmark(namespace):
    if namespace=='BPO':
        ancestor_path = '/home/nzhou/git/CAFAAssess/precrec/gene_ontology_edit.obo_ancestors_bpo.txt'
        benchmark_path = '/home/nzhou/git/CAFAAssess/precrec/leafonly_BPO.txt'
    elif namespace=='MFO':
        ancestor_path = '/home/nzhou/git/CAFAAssess/precrec/gene_ontology_edit.obo_ancestors_mfo.txt'
        benchmark_path = '/home/nzhou/git/CAFAAssess/precrec/leafonly_MFO.txt'
    else:
        ancestor_path = '/home/nzhou/git/CAFAAssess/precrec/gene_ontology_edit.obo_ancestors_cco.txt'
        benchmark_path = '/home/nzhou/git/CAFAAssess/precrec/leafonly_CCO.txt'
                                                                             
    bench = benchmark(ancestor_path,benchmark_path)
    bench.propagate()
    return bench

      
'''            
The top ten teams in BPO:  
BPO:          
Paccanaro Lab:117
Tian Lab:85
MS-kNN:129
Jones-UCL:114,120,125
Orengo-FunFams: 127
Gough Lab:110,111,112
Argot2:94
PULP: 83
PFPDB:75
INGA-Tosatto:119 (no non-HPO human)

MFO:
MS-knn: 129
*EVEX:72
tianlab: 85
PFPDB:75
Orengo-FunFams: 127
*Go2Proto:101
*SIFTER 2.4: 102
INGA-Tosatto:119 (no non-HPO human)

Jones-UCL: 
Argot2: 94

CCO:
EVEX:72
Tian Lab:85
Gough Lab:
MS-KNN:129
JONES-UCL
*CONS:115
*Rost Lab:97
PULP:83
*IASL:86
ORENGO-FUNFAMS:127

Naive
BLAST:        
'''  

def prediction_ontology_split_write(pred_path, obo_path):
    """
    Separate the prediction file into the different ontologies
    pred_path should be a handle!!!!!!!!!!
    This should work
    Last edited by Ashley 05/11/2016
    
    This is edited to change the output file directory! 05/23/2016
    This function is different from the one in precRec.py!!
    """
    
    all_pred = GOPred()
    #pred_path = open('/home/nzhou/git/CAFAAssess/precrec/M1HS.74.Homo_sapiens.txt')
    all_pred.read(pred_path)
    #obo_path = '/home/nzhou/git/Ontology/go-basic.obo'
    go_graph = OboIO.OboReader(open(obo_path)).read()
    mfo_out = open("%s_MFO.txt" % pred_path.name.split('/')[-1].split('.')[-2],"w")
    bpo_out = open("%s_BPO.txt" % pred_path.name.split('/')[-1].split('.')[-2],"w")
    cco_out = open("%s_CCO.txt" % pred_path.name.split('/')[-1].split('.')[-2],"w")
    for protein in all_pred.data.items():
        for u in protein[1]:
            if go_graph.get_namespace(u['term']) == 'molecular_function': 
                mfo_out.write("%s\t%s\t%.2f\n" % (protein[0], u['term'], u['confidence']))
            elif go_graph.get_namespace(u['term']) == 'biological_process': 
                bpo_out.write("%s\t%s\t%.2f\n" % (protein[0], u['term'], u['confidence']))
            elif go_graph.get_namespace(u['term']) == 'cellular_component': 
                cco_out.write("%s\t%s\t%.2f\n" % (protein[0], u['term'], u['confidence']))
            else:
                raise ValueError ("Term %s not found in any ontology" % u['term'])
    mfo_out.close()
    bpo_out.close()
    cco_out.close()



from os import walk
def getFileNames(pred_folder):
    filenames = []
    newfilenames = []
    for (dirpath,dirnames,filename) in walk(pred_folder):
        filenames.extend(filename)
        break
    return filenames 

def pred_split(predfile,submission_folder):
    #sys.path.append(submission_folder)
    if predfile.split('.')[-1]=="zip" :
        obj = zipfile.ZipFile(submission_folder+predfile,'r')
        #obj = zipfile.ZipFile('/home/nzhou/old computer/Documents/CAFA2/CAFA2_submissions/83/M2HS.83.bonneauLab_2_9606.txt.zip')
        for filename in obj.namelist():
            if 'hpo' not in filename and filename.split('.')[-1]=='txt':
                 handle = obj.open(filename,'r')
                 print filename
                 prediction_ontology_split_write(handle,obo_path)
                 handle.close()
        obj.close()
    elif predfile.split('.')[-1]=='tgz' or ''.join(predfile.split('.')[-2:])=='targz':
        obj = tarfile.open(submission_folder+predfile,'r')
        for filename in obj.getnames():
            if 'hpo' not in filename and filename.split('.')[-1]=='txt':
                handle = obj.extractfile(filename)
                prediction_ontology_split_write(handle,obo_path)
                handle.close()
        obj.close() 
    elif predfile.split('.')[-1]=="gz":
        handle = gzip.open(submission_folder+predfile,'r')
        prediction_ontology_split_write(handle,obo_path)
        handle.close()
    elif predfile.split('.')[-1]=="txt":
        handle = open(submission_folder+predfile,'r')
        prediction_ontology_split_write(handle,obo_path)
        handle.close()
       
def files_split(teamNumber,human):
    '''
    human is boolean
    if true, only work on human prediction file
    '''
    submission_folder = '/home/nzhou/old computer/Documents/CAFA2/CAFA2_submissions/'+str(teamNumber)+'/'
    filelist = getFileNames(submission_folder)
    #os.chdir(submission_folder)
    pred_folder = '/home/nzhou/git/CAFAAssess/confidence/'+str(teamNumber)+"/" 
    if not os.path.exists(pred_folder):
        os.mkdir(pred_folder)
    os.chdir(pred_folder)
    for i in filelist:
        if human:
            if '9606' in i or 'sapien' in i:
                print i
                pred_split(i,submission_folder)
                continue
            else:
                continue
        else:
            pred_split(i,submission_folder)

def get_namespace_index(namespace):
    num = None
    if namespace=='BPO' or namespace=='bpo':
        num = 0
    elif namespace=='MFO' or namespace=='mfo':
        num = 1
    elif namespace=='CCO' or namespace=='cco':
        num =2
    else:
        raise ValueError("name space name not found, check prediction files")
        print namespace
    return num
            

'''
main
'''
#


if __name__ == '__main__':
    obo_path = '/home/nzhou/git/CAFAAssess/precrec/gene_ontology_edit.obo.2014-06-01'
    #teamNums = [117,85,129,127,94,83,75]
    #teamNums = [85,129,127,94,83,75]
    teamNums = [115,97,86]
    #teamNums = [72,101,102,115,97,86]
    #These teams were split to 
    bench = [read_benchmark('BPO'), read_benchmark('MFO'), read_benchmark('CCO')]
    #For team 115, prediction file is bad, with "GO:GO:", So the orginal file was untarred in the submission folder
    for num in teamNums:   
        if num!=115:
            files_split(num,True)
        os.chdir('/home/nzhou/git/CAFAAssess/confidence/'+str(num)+'/')
        files = getFileNames(os.getcwd())
        for f in files:
            if 'confdata' not in f and f.split('.')[-1]=='txt' and os.path.getsize(f)!=0:
                print f
                pred = GOPred()
                pred.read(open(f))
                namespace = f.split('.')[-2][-3:]
                b = bench[get_namespace_index(namespace)]
                pr = PrecREC(b,pred)
                print os.getcwd() 
                pr.printConfidence(os.getcwd()+'/'+f.split('.')[0]+'_confdata.txt')
'''
    baseline = ["BLAST","Naive"]
    for name in baseline:
        os.chdir('/home/nzhou/git/CAFAAssess/confidence/'+name+'/')
        files = getFileNames(os.getcwd())
        for f in files:
            if 'hpo' in f:
                continue
            if 'confdata' not in f and f.split('.')[-1]=='txt':
                pred = GOPred()
                pred.read(open(f))
                namespace = f.split('.')[-2][-3:]
                b = bench[get_namespace_index(namespace)]
                pr = PrecREC(b,pred)
                pr.printConfidence(os.getcwd()+'/'+f.split('.')[0]+'_confdata.txt')
'''  


