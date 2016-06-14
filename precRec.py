# -*- coding: utf-8 -*-
#!/usr/bin/env python3
"""
new CAFA precision recall assessment

Created on Thu Apr 14 17:14:30 2016
dependency: - 
@author: Ashley Zhou
"""



import os
import sys
import numpy
import matplotlib.pyplot as plt
from collections import defaultdict
#os.chdir('/home/nzhou/git')
#from Ontology.IO import OboIO
print os.getcwd()
sys.path.append('/home/nzhou/git')
from CAFAAssess.precrec.GOPred import GOPred


class benchmark:
    def __init__(self,ancestor_path,benchmark_path):
        '''
        Here benchmark_path should be ontology specific
        ancestor_path should also be ontology specific
        '''
        #Key: protein
        #Value: set of benchmark leaf terms
        self.ancestors = defaultdict(set)
        # Read GO ancestors file generated with go_ontology_ancestors_split_write()
        # File format: 
        # go_term <tab> ancestor_1,ancestor_2,..,ancestor_n
        with open(ancestor_path) as ancestors_input:
            for inline in ancestors_input:
                inrec = inline.strip().split('\t')
                term = inrec[0]
                if len(inrec) == 1:
                    self.ancestors[term] = set({})
                else:
                    term_ancestors = inrec[1]
                    self.ancestors[term] = set(term_ancestors.split(','))
                
        self.true_base_terms = defaultdict(set)
        with open(benchmark_path) as benchmark_input:
            for inline in benchmark_input:
                protein, term = inline.strip().split('\t')
                self.true_base_terms[protein].add(term)
                
    def propagate(self):
        #Key: protein
        #Value: set of benchmark propagated terms
        self.true_terms = defaultdict(set)
        for protein in self.true_base_terms:
            for term in self.true_base_terms[protein]:
                try:
                    ancestors = self.ancestors[term]
                except KeyError:
                    sys.stderr.write("%s not found \n" % term) 
                self.true_terms[protein].add(term)
                #A: what does |= do?
                self.true_terms[protein] |= ancestors
        
    

class PrecREC:
    '''
    New code by Ashley
    updated: 05/24/2016
    A class for doing precision recall calculations
    Not necssarily curve computation
    read in ONTOLOGY-SPECIFIC prediction file
    #Do we need ONTOLOGY-SPECIFIC ontology file??
    '''
    def __init__(self, benchmark, GOPred):
        '''
        constructor
        benchmark is an instance of the benchmark class
        countb is the number of predicted proteins in this file that are in the benchmark file
        counta is the number of proteins with at least one term above threshold
        '''
        self.ancestors = benchmark.ancestors
        self.true_terms = benchmark.true_terms
        self.obsolete = set()
        self.counta = defaultdict()
        self.countb = 0
        
        #predicted_base_terms is the same as GOPred.data
        #No need creating another dictionary
        
        #Now propogate the predicted terms
        #key:protein
        #value: list of dictionaries
        #key: GO term
        #value: tuple(confidence, True/False) whether in true terms or not
        #take the largest confidence
        #Take care of obsolete terms as well
        
        #
        self.predicted = defaultdict(defaultdict)
        for prot in GOPred.data:
            if benchmark.true_terms[prot]!=set():
                '''
                The protein is in the benchmark file
                i.e. gained experimental annotation
                '''
                self.countb += 1
                for tc in GOPred.data[prot]:
                    try:
                        ancterms = self.ancestors[tc['term']]
                    except KeyError:
                        #sys.stderr.write("%s not found\n" % tc['term'])
                        self.obsolete.add(tc['term'])
                        continue
                    if ancterms == set() and tc['term'] not in (u'GO:0003874', u'GO:0008150' , u'GO:0005575'):
                        #sys.stderr.write("no ancestors found for %s\n" % tc['term'])
                        self.obsolete.add(tc['term'])
                        continue
                    
                    if tc['term'] in self.predicted[prot]:
                        #This term has already been added
                        #maybe as an ancestor of other terms
                        #update confidence with the maximum one
                        #propagate and update all ancestor confidence
                        self.__update_confidence__(prot,tc)
                    else:
                        #add this term to self.predicted
                        #add confidence, and compare with self.true_terms
                        #if term in self.true_terms, True
                        #else False
                        #No matter true or false, propagate
                        self.predicted[prot][tc['term']]=self.__compare__(prot,tc)
                        for ancterm in ancterms:
                            newtc = {'term':ancterm,'confidence':tc['confidence']}
                            if ancterm in self.predicted[prot]:
                                self.__update_confidence__(prot,newtc)
                            else:
                                self.predicted[prot][ancterm]=self.__compare__(prot,newtc)
            else:
                '''
                this protein is not in the benchmark file
                '''
                self.predicted[prot] = None
                


    def __update_confidence__(self,prot,tc):
        '''
        tc is a dictionary with {'confidence':0.57,'term':'GO:006644'}
        prot is a protein
        This function compares the confidence value in tc, and if it's larger than
        the confidence that's been added for term in self.predicted
        we update that confidence
        It also updates all propagated terms of tc
        '''
        if tc['confidence']>self.predicted[prot][tc['term']][0]:
            self.predicted[prot][tc['term']][0]=tc['confidence']
            for ancterm in self.ancestors[tc['term']]:
                if tc['confidence']>self.predicted[prot][ancterm][0]:
                    self.predicted[prot][ancterm][0]=tc['confidence']
                    
                    
    def __compare__(self,prot,tc):
        '''
        tc is a dictionary with {'confidence':0.57,'term':'GO:006644'}
        prot is a protein
        This function compares if tc['term'] is in self.true_terms
        returns a list ["confidence","True/False"]
        '''
        if tc['term'] in self.true_terms[prot]:
            return [tc['confidence'],True]
        else:
            return [tc['confidence'],False]
            
        
    def getObsolete(self):
        '''
        return all obsolete terms used by the prediction team
        '''
        return(self.obsolete)

    def term_precision_recall(self,threshold,protein):
        TP = 0.0
        count = 0
        #count is to count how many terms are above the threshold
        if self.predicted[protein] is not None:
            for term in self.predicted[protein]:
                if self.predicted[protein][term][0]>=threshold:
                    #greater but not greater or equal
                    count+=1
                    if self.predicted[protein][term][1] :
                        TP+=1
            try:
                precision = TP/count
            except ZeroDivisionError:
                precision=None
            recall = TP/len(self.true_terms[protein])
            #recall should not have zerodivision problem
            #since if self.predicted[protein] is not None
            #This protein is in the benchmark file
            #i.e. gained experimental annotation
            #len(self.true_terms[protein]) should not be 0
            return (precision, recall)
        else:
            return (None,None)
                        
    def precision_recall(self,threshold):
        '''
        this calculates the overall precision recall of the team, given a threshold,
        For one prediction file, i.e. for one species and one model!!!!
        05/24/2016: countb has been calculated over and over again for each threshold
        '''
        prec = float(0)
        self.counta[threshold] = 0
        rec = float(0)
        for prot in self.predicted:
                a,b = self.term_precision_recall(threshold,prot)
                if a is not None:
                    prec +=a
                    self.counta[threshold]+=1
                if b is not None:
                    rec +=b
        return (prec/self.counta[threshold], rec/self.countb)
    
    def getNumProteins(self,threshold):
         
        '''
        run precision_recall first
        '''         
        print('number of benchmark proteins: %s\n'% len(self.true_terms))
        print('number of proteins predicted: %s\n'% len(self.predicted))
        #Those with not-None recall: 
        #(predicted protein that are in benchmark, only one species per prediction file!! )
        print ('number of predicted proteins that are in the benchmark file: %s\n' % self.countb)
        #Those with not-None precision:
        try:
            print('number of proteins with at least one term above threshold: %s\n' % self.counta[threshold] )
        except KeyError:
            sys.stderr.write("Run precision_recall(%s) first\n" % str(threshold))
            
    def Fmax_output(self,interval):
        '''
        returns the fmax value AND outputs the precision-recall values for each threshold
        This computes precision and recall for every threshold
        interval is number of threshold values between 0 and 1:
        e.g. interval = 9, thresholds: 0.   ,  0.125,  0.25 ,  0.375,  0.5  ,  0.625,  0.75 ,  0.875,  1. 
        '''
        fmax = 0
        pre = []
        rec = []
        for thres in numpy.linspace(0.01,0.99,interval):
            a,b = self.precision_recall(thres)
            pre.append(a)
            rec.append(b)
            f = 2*a*b/(a+b)
            if f>=fmax:
                fmax = f
        return (pre,rec,fmax)
    
    def printConfidence(self,output_path):
        '''
        print confidence and True/False to a file
        to be read by R
        '''
        protindex = 0
        out = open(output_path,'w')
        for prot in self.predicted:
            protindex += 1
            if self.predicted[prot] is not None:
                for term in self.predicted[prot]:
                    out.write("%s\t%s\t%s\n" % (str(protindex), self.predicted[prot][term][0],self.predicted[prot][term][1]))
        out.close()

        
        




if __name__=='__main__':
    '''
    ancestor_path = '/home/nzhou/git/CAFAAssess/precrec/gene_ontology_edit.obo_ancestors_bpo.txt'  
    #This benchmark is the full BPO benchmark
    #benchmark_path = '/home/nzhou/git/CAFAAssess/precrec/leafonly_BPO.txt'
    #bench = benchmark(ancestor_path,benchmark_path)
    #bench.propagate()
    #This benchmark is the human only BPO benchmark
    hs = benchmark(ancestor_path,'/home/nzhou/git/CAFAAssess/precrec/leafonly_BPO_9606.txt')
    hs.propagate()
    
    all_pred = GOPred()
    pred_path ='/home/nzhou/git/CAFAAssess/confidence/117/PaccanaroLab_1_9606_BPO.txt'
    all_pred.read(open(pred_path))
    #b = PrecREC(bench,all_pred)
    c = PrecREC(hs,all_pred)
    #c.printConfidence('/home/nzhou/git/CAFAAssess/confidence/117/newconfdata.txt')
    #print b.precision_recall(0.1)
    print c.precision_recall(0.1)
    #print b.getNumProteins(0.1)
    #print c.getNumProteins(0.1)
    fm = c.Fmax_output(9)
    plt.plot(fm[0],fm[1],'ro')
    plt.axis([0,1,0,1])
    plt.show()
    '''
    '''
    get confidence files for all human predictions for all top ten teams (model 1?)
    in confidence.py
    '''
    '''
    Get PR curve and Fmax for all top-performing teams in three ontologies
    '''
    bpo = ('/117/PaccanaroLab_1_9606_BPO.txt','/85/tianlab_1_9606_BPO.txt','/129/tu_1_9606_BPO.txt','/127/txt_BPO.txt','/94/txt_BPO.txt','/83/bonneauLab_1_9606_BPO.txt','/75/Homo_sapiens_BPO.txt','/Naive/BN4S_bpo.txt','/BLAST/BB4S_bpo.txt')
    mfo = ('/129/tu_1_9606_MFO.txt','/72/EVEX_1_9606_MFO.txt','/85/tianlab_1_9606_MFO.txt','/75/Homo_sapiens_MFO.txt','/127/txt_MFO.txt','/101/GO2Proto_1_9606_0_MFO.txt','/102/brennersifter_1_9606_0_MFO.txt','/94/txt_MFO.txt','/Naive/BN4S_mfo.txt','/BLAST/BB4S_mfo.txt')
    cco = ('/72/EVEX_1_9606_CCO.txt','/85/tianlab_1_9606_CCO.txt','/129/tu_1_9606_CCO.txt','/115/Homo_sapiens_CCO.txt','/83/bonneauLab_1_9606_CCO.txt','/86/IASLAcademiaSinica_1_9606_CCO.txt','/127/txt_CCO.txt','/Naive/BN4S_cco.txt','/BLAST/BB4S_cco.txt')
    
    #before running this function, make sure that a benchmark instance has been created and propagated.
    def getplot(filename):
        fullpath = '/home/nzhou/git/CAFAAssess/confidence'+filename
        all_pred = GOPred()
        all_pred.read(open(fullpath))
        c = PrecREC(hs,all_pred)
        fm = c.Fmax_output(99)
        
        plt.plot(fm[1],fm[0])
        plt.axis([0,1,0,1])
        plt.yticks(numpy.arange(0,1,0.1))
        plt.xlabel('Recall')
        plt.ylabel('Precision')
        plt.title(fullpath.split('/')[-1].split('.')[0])
        plt.savefig(fullpath.split('.')[0]+'.png')
        plt.close()
        return fm[2]
        
    #BPO
    ancestor_path = '/home/nzhou/git/CAFAAssess/precrec/gene_ontology_edit.obo_ancestors_bpo.txt'  
    benchmark_path =  '/home/nzhou/git/CAFAAssess/precrec/leafonly_BPO.txt'
    hs = benchmark(ancestor_path, benchmark_path)
    hs.propagate()
    bpofm = defaultdict()
    for i in bpo:
        bpofm[i] = getplot(i)
        
    #MFO
    ancestor_path = '/home/nzhou/git/CAFAAssess/precrec/gene_ontology_edit.obo_ancestors_mfo.txt'  
    benchmark_path =  '/home/nzhou/git/CAFAAssess/precrec/leafonly_MFO.txt'
    hs = benchmark(ancestor_path, benchmark_path)
    hs.propagate()
    mfofm = defaultdict()
    for j in mfo:
        mfofm[j]=getplot(j)
    #CCO
    ancestor_path = '/home/nzhou/git/CAFAAssess/precrec/gene_ontology_edit.obo_ancestors_cco.txt'  
    benchmark_path =  '/home/nzhou/git/CAFAAssess/precrec/leafonly_CCO.txt'
    hs = benchmark(ancestor_path, benchmark_path)
    hs.propagate()
    ccofm = defaultdict()
    for k in cco:
        ccofm[k]=getplot(k)
        

                    
    
    
