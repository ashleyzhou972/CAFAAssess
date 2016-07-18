
# A command-line interface for assessing a single CAFA prediction file using precision-recall.

All GO files, benchmark files and test prediction files are in /precrec/. 

Default GO release is from 06/01/2014.

To use:

### 1. Download both CAFAAssess and Ontology repositories

### 2. cd to the nearest base directory that contains both CAFAAssess and Ontology

### 3. Type `python CAFAAssess/precrec_main.py -h` for the usage of this module

## Example:

`python CAFAAssess/precrec_main.py BPO 117 9606 ./CAFAAssess/TEST_1_9606.txt ./CAFAAssess/testplot.png`


Future work including another module that processes zipped submission files and combining species
