#!/usr/bin/env python3

#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>

import re
import sys
from collections import defaultdict


pr_field = re.compile("^PR=[0,1]\.[0-9][0-9];$")
rc_field = re.compile("^RC=[0,1]\.[0-9][0-9]$")
# Fix to add HP ontology 2014-1-9
go_field = re.compile("^(GO|HP):[0-9]{5,7}$")
# Fix to add EFI targets 2014-1-9
target_field = re.compile("^(T|EFI)[0-9]{5,20}$")
confidence_field = re.compile("^[0,1]\.[0-9][0-9]$")

# Legal states: the CAFA prediction records fields, and their order. KEYWORDS and ACCURACY are
# optional
legal_states1 = ["author","model","keywords","accuracy","go_prediction","end"]
legal_states2 = ["author","model","keywords","go_prediction","end"]
legal_states3 = ["author","model","go_prediction","end"]

legal_keywords = [
"sequence alignment", "sequence-profile alignment", "profile-profile alignment", "phylogeny",
"sequence properties",
"physicochemical properties", "predicted properties", "protein interactions", "gene expression",
"mass spectrometry",
"genetic interactions", "protein structure", "literature", "genomic context", "synteny", 
"structure alignment",
"comparative model", "predicted protein structure", "de novo prediction", "machine learning", 
"genome environment", 
"operon", "ortholog", "paralog", "homolog", "hidden markov model", "clinical data", "genetic data", 
"natural language processing", "other functional information"
]
class GOPred:
    """
    A class for reading and storing CAFA GO predictions
    self.data is a dictionary
       key: protein ID
       value: [{'term':go_term_1, 'threshold': threshold_1},...,{'term':go_term_n, 'threshold': threshold_n}]
    """
    def __init__(self):
        self.author = None
        self.model = None
        self.keywords = []
        self.data = defaultdict(list)

    def _author_check(self,inrec):
        correct = True
        errmsg = None
        fields = [i.strip() for i in inrec.split()]
        if len(fields) != 2:
            correct = False
            errmsg = "AUTHOR: invalid number of fields. Should be 2"
        elif fields[0] != "AUTHOR":
            correct = False
            errmsg = "AUTHOR: First field should be AUTHOR"
        else:
            self.author = fields[1]
        return correct, errmsg

    def _model_check(self,inrec):
        correct = True
        errmsg = None
        fields = [i.strip() for i in inrec.split()]
        if len(fields) != 2:
            correct = False
            errmsg = "MODEL: invalid number of fields. Should be 2"
        elif fields[0] != "MODEL":
            correct = False
            errmsg = "MODEL: First field should be MODEL"
        elif len(fields[1]) != 1 or not fields[1].isdigit():
            correct = False
            errmsg = "MODEL: second field should be single digit."
        else:
            self.model = int(fields[1])
        return correct, errmsg

    def _keywords_check(self,inrec):
        correct = True
        errmsg = None
        if inrec[:8] != "KEYWORDS":
            correct = False
            errmsg = "KEYWORDS: first field should be KEYWORDS"
        else:
            keywords = [i.strip().lower() for i in inrec[8:].split(",")]
            for keyword in keywords:
                #stupid full stop 
                if keyword[-1] == ".":
                    keyword = keyword[:-1]
                if keyword not in legal_keywords:
                    correct = False
                    errmsg = "KEYWORDS: illegal keyword %s" % keyword
                    break
                else:
                    self.keywords.append(keyword)
        return correct, errmsg

    def _accuracy_check(self,inrec):
        correct = True
        errmsg = None
        fields = [i.strip() for i in inrec.split()]
        if len(fields) != 4:
            correct = False
            errmsg = "ACCURACY: error in number of fields. Should be 4"
        elif fields[0] != "ACCURACY":
            correct = False
            errmsg = "ACCURACY: first field should be 'ACCURACY'"
        elif not fields[1].isdigit() or len(fields[1]) != 1:
            correct = False
            errmsg = "ACCURACY: second field should be a single digit"
        elif not pr_field.match(fields[2]):
            correct = False
            errmsg = "ACCURACY: error in PR field"
        elif not rc_field.match(fields[3]):
            correct = False
            errmsg = "ACCURACY: error in RC field"
        return correct, errmsg


    def _go_prediction_check(self,inrec):
        correct = True
        errmsg = None
        fields = [i.strip() for i in inrec.split()]
        if len(fields) != 3:
            correct = False
            errmsg = "GO prediction: wrong number of fields. Should be 3"
        elif not target_field.match(fields[0]):
            correct = False
            errmsg = "GO prediction: error in first (Target ID) field"
        elif not go_field.match(fields[1]):
            correct = False
            errmsg = "GO prediction: error in second (GO ID) field"
        elif not confidence_field.match(fields[2]):
            correct = False
            errmsg = "GO prediction: error in third (confidence) field"
        elif float(fields[2]) > 1.0:
            correct = False
            errmsg = "GO prediction: error in third (confidence) field. Cannot be > 1.0"
        else:
            self.data[fields[0]].append({'term': fields[1], 'confidence': float(fields[2])})
        return correct, errmsg

    def _end_check(self,inrec):
        correct = True
        errmsg = None
        fields = [i.strip() for i in inrec.split()]
        if len(fields) != 1:
            correct = False
            errmsg = "END: wrong number of fields. Should be 1"
        elif fields[0] != "END":
            correct = False
            errmsg = "END: record should include the word END only"
        return correct, errmsg

    def _handle_error(self,correct, errmsg, inrec):
        if not correct:
            print (inrec)
            raise ValueError(errmsg)


    def read(self, inpath):
        visited_states = []
        s_token = 0
        n_accuracy = 0
        first_prediction = True
        first_accuracy = True
        first_keywords = True
        n_models = 0

        with open(inpath) as inpred:
            for inline in inpred: 
                inrec = [i.strip() for i in inline.split()]
                field1 = inrec[0]
                # Check which field type (state) we are in
                if field1 == "AUTHOR":
                    state = "author"
                elif field1 == "MODEL":
                    state = "model"
                elif field1 == "KEYWORDS":
                    state = "keywords"
                elif field1 == "ACCURACY":
                    state = "accuracy"
                elif field1 == "END":
                    state = "end"
                else: #default to prediction state
                    state = "go_prediction"
                # Check for errors according to state
                if state == "author":
                    correct,errmsg = self._author_check(inline)
                    self._handle_error(correct, errmsg,inline)
                    visited_states.append(state)
                elif state == "model":
                    n_models += 1
                    n_accuracy = 0
                    if n_models > 3:
                        raise ValueError("Too many models. Only up to 3 allowed")
                    correct,errmsg = self._model_check(inline)
                    self._handle_error(correct, errmsg,inline)
                    if n_models == 1:
                        visited_states.append(state)
                elif state == "keywords":
                    if first_keywords:
                        visited_states.append(state)
                        first_keywords = False
                    correct, errmsg = self._keywords_check(inline)
                    self._handle_error(correct, errmsg,inline)
                elif state == "accuracy":
                    if first_accuracy:
                        visited_states.append(state)
                        first_accuracy = False
                    n_accuracy += 1
                    if n_accuracy > 3:
                        self._handle_error(False, "ACCURACY: too many ACCURACY records")
                    else:
                        correct, errmsg = self._accuracy_check(inline)
                elif state == "go_prediction":
                    correct, errmsg = self._go_prediction_check(inline)
                    self._handle_error(correct, errmsg,inline)
                    if first_prediction:
                        visited_states.append(state)
                        first_prediction = False
                elif state == "end":
                    correct, errmsg = self._end_check(inline)
                    self._handle_error(correct, errmsg,inline)
                    visited_states.append(state)
            # End file forloop
            if (visited_states != legal_states1 and
                visited_states != legal_states2 and
                visited_states != legal_states3):
                print (visited_states)
                print ("file not formatted according to CAFA specs")
                print ("Check whether all these record types are in your file")
                print ("Check whether all these record types are in your file in the correct order")
                print ("AUTHOR, MODEL, KEYWORDS, ACCURACY (optional), predictions, END")
                raise ValueError

class PrecRec:
    """
    Precision-recall calculations for CAFA
    """
    def __init__(self):
        self.ancestors = {}

        # Key: protein; Value: set of true terms from CAFA benchmark
        self._true_base_terms = defaultdict(set)

        # Key: protein; value: set of terms from CAFA benchmark, and their ancestors
        self.true_terms = defaultdict(set)

    def read_ancestors(self,ancestors_path):
        # Call this method first, populates self.ancestors
        # Read GO ancestors file generated with obo2ancestors
        # File format: 
        # go_term <tab> ancestor_1,ancestor_2,..,ancestor_n
        with open(ancestors_path) as ancestors_input:
            for inline in ancestors_input:
                inrec = inline.strip().split('\t')
                term = inrec[0]
                if len(inrec) == 1:
                    self.ancestors[term] = set({})
                else:
                    term_ancestors = inrec[1]
                    self.ancestors[term] = set(term_ancestors.split(','))

    def read_benchmark(self, benchmark_path):
        # Read the benchmark file that contains the true annotations
        # 
        with open(benchmark_path) as benchmark_input:
            for inline in benchmark_input:
                protein, term = inline.strip().split('\t')
                self._true_base_terms[protein].add(term)

    def propagate_true_terms(self):
        for protein in self._true_base_terms:
            for term in self._true_base_terms[protein]:
                try:
                    ancestors = self.get_ancestors(term)
                except KeyError:
                    print("not found %s" % term) 
                self.true_terms[protein].add(term) # Add the base term
                self.true_terms[protein] |= ancestors # Add all ancestors

    def get_ancestors(self,term):
        return self.ancestors[term]
    
    def term_precision_recall(self, protein, term):
        # Precision-recall for a single term for a single protein
        # term: single GO term to be propagated
        # true_terms: set of true terms

        pred_terms = self.get_ancestors(term)
        pred_terms.add(term)
        prot_true_terms = self.true_terms[protein]
        precision = len(pred_terms & prot_true_terms) / len(pred_terms)
        recall = len(pred_terms & prot_true_terms) / len(prot_true_terms)
        return (precision, recall)

def precision_recall(prediction, benchmark_path, ancestors_path):
    # accepts a GOPred instantiation
    # does precision / recall calculation
    # TODO: take care of threshold!!
    prec_rec_vector = []
    prec_rec = PrecRec()
    prec_rec.read_ancestors(ancestors_path)
    prec_rec.read_benchmark(benchmark_path)
    for threshold in [i*0.01 for i in range(1,101)]:
        for protein in prediction.data:
            for term_thresh in prediction.data[protein]:
                term = term_thresh['term']
                thresh_level = term_thresh['threshold']
                if threhsold >= thresh_level:
                    precision, recall = prec_rec.term_precision_recall(protein, term)
                    prec_rec_vector.append((precision, recall))
    return prec_rec_vector


if __name__ == '__main__':
    my_gopred = GOPred()
    my_gopred.read(sys.argv[1]) # Read a prediction
    for i in list(my_gopred.data.items())[:10]:
        print(i)

    
