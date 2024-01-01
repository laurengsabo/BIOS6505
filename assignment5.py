"""
Assignment 5: Automatic Primer Identification

In this assignment, you will practice using python's subprocess library, the pyfasta library, and the primer3
command-line-interface to automatically identify primers. This template provides a potential route to a completed
script, but the implementation details are entirely up to you. The checker script, as you will see, does little more
than confirm that your script prints the correct primers to the console for a handful of argument combinations.
Author: Tucker J Lancaster
Modified by: Lauren Sabo and Karthik Krishnan 
"""

import argparse, pdb
import subprocess
from pyfasta import Fasta

"""
If you have not already, install pyfasta and primer3 into the current environment using the following commands:
conda install -c bioconda primer3
conda install -c bioconda pyfasta
"""

"""
Create your parser and arguments. Your script should expect three positional arguments: First, the name of a
fasta file, second, the chromosome you want to amplify, and third, the position you want amplified. Enforce that
the first two arguments are of type str, and the third is of type int.
"""

parser = argparse.ArgumentParser(description="Automatically identify primers.")
parser.add_argument("fasta_file", help="Name of the fasta file", type=str)
parser.add_argument("chromosome", help="Chromosome to amplify", type=str)
parser.add_argument("position", help="Position to amplify", type=int)
parser.add_argument("--product_length", type=int, default=700, help="Product length (default: 700)")
args = parser.parse_args()

"""
Open the fasta file and read in the entirety of the sequence. You can use whatever data structure you prefer.
"""

fasta = Fasta(args.fasta_file)

"""
Identify the sequence 500 bp upstream and downstream of the requested position and store it as a string.
"""

start = args.position - 500
end = args.position + 500
req_seq = str(fasta[args.chromosome][start:end])
print(req_seq)

"""
Add braces to the DNA sequence 100 bp upstream and downstream of the requested position. Alternatively, you can use 
an argument in the input file to specify this
"""

#put in as an argument

"""
Create an input file containing the DNA sequence and specify a product length of 600 - 800 base pairs
"""
 
with open('primer.txt', 'w') as new_file:
    new_file.write("SEQUENCE_TEMPLATE=" + req_seq + "\n")
    new_file.write("PRIMER_ANNEALING_TEMP=52.0" + "\n")
    new_file.write("PRIMER_DMSO_CONC=0.0" + "\n")
    new_file.write("PRIMER_DMSO_FACTOR=0.6" + "\n")
    new_file.write("PRIMER_DNA_CONC=50.0" + "\n")
    new_file.write("PRIMER_DNTP_CONC=0.6" + "\n")
    new_file.write("PRIMER_EXPLAIN_FLAG=1" + "\n")
    new_file.write("PRIMER_FIRST_BASE_INDEX=1" + "\n")
    new_file.write("PRIMER_FORMAMIDE_CONC=0.0" + "\n")
    new_file.write("PRIMER_GC_CLAMP=0" + "\n")
    new_file.write("PRIMER_INSIDE_PENALTY=-1.0" + "\n")
    new_file.write("PRIMER_INTERNAL_DMSO_CONC=0.0" + "\n")
    new_file.write("PRIMER_INTERNAL_DMSO_FACTOR=0.6" + "\n")
    new_file.write("PRIMER_INTERNAL_DNA_CONC=50.0" + "\n")
    new_file.write("PRIMER_INTERNAL_DNTP_CONC=0.0" + "\n")
    new_file.write("PRIMER_INTERNAL_FORMAMIDE_CONC=0.0" + "\n")
    new_file.write("PRIMER_INTERNAL_MAX_BOUND=110.0" + "\n")
    new_file.write("PRIMER_INTERNAL_MAX_GC=80.0" + "\n")
    new_file.write("PRIMER_INTERNAL_MAX_HAIRPIN_TH=47.00" + "\n")
    new_file.write("PRIMER_INTERNAL_MAX_LIBRARY_MISHYB=12.00" + "\n")
    new_file.write("PRIMER_INTERNAL_MAX_NS_ACCEPTED=0" + "\n")
    new_file.write("PRIMER_INTERNAL_MAX_POLY_X=5" + "\n")
    new_file.write("PRIMER_INTERNAL_MAX_SELF_ANY=12.00" + "\n")
    new_file.write("PRIMER_INTERNAL_MAX_SELF_ANY_TH=47.00" + "\n")
    new_file.write("PRIMER_INTERNAL_MAX_SELF_END=12.00" + "\n")
    new_file.write("PRIMER_INTERNAL_MAX_SELF_END_TH=47.00" + "\n")
    new_file.write("PRIMER_INTERNAL_MAX_SIZE=27" + "\n")
    new_file.write("PRIMER_INTERNAL_MAX_TM=63.0" + "\n")
    new_file.write("PRIMER_INTERNAL_MIN_3_PRIME_OVERLAP_OF_JUNCTION=4" + "\n")
    new_file.write("PRIMER_INTERNAL_MIN_5_PRIME_OVERLAP_OF_JUNCTION=7" + "\n")
    new_file.write("PRIMER_INTERNAL_MIN_BOUND=-10.0" + "\n")
    new_file.write("PRIMER_INTERNAL_MIN_GC=20.0" + "\n")
    new_file.write("PRIMER_INTERNAL_MIN_QUALITY=0" + "\n")
    new_file.write("PRIMER_INTERNAL_MIN_SIZE=18" + "\n")
    new_file.write("PRIMER_INTERNAL_MIN_THREE_PRIME_DISTANCE=-1" + "\n")
    new_file.write("PRIMER_INTERNAL_MIN_TM=57.0" + "\n")
    new_file.write("PRIMER_INTERNAL_OPT_BOUND=97.0" + "\n")
    new_file.write("PRIMER_INTERNAL_OPT_GC_PERCENT=50.0" + "\n")
    new_file.write("PRIMER_INTERNAL_OPT_SIZE=20" + "\n")
    new_file.write("PRIMER_INTERNAL_OPT_TM=60.0" + "\n")
    new_file.write("PRIMER_INTERNAL_SALT_DIVALENT=0.0" + "\n")
    new_file.write("PRIMER_INTERNAL_SALT_MONOVALENT=50.0" + "\n")
    new_file.write("PRIMER_INTERNAL_WT_BOUND_GT=0.0" + "\n")
    new_file.write("PRIMER_INTERNAL_WT_BOUND_LT=0.0" + "\n")
    new_file.write("PRIMER_INTERNAL_WT_END_QUAL=0.0" + "\n")
    new_file.write("PRIMER_INTERNAL_WT_GC_PERCENT_GT=0.0" + "\n")
    new_file.write("PRIMER_INTERNAL_WT_GC_PERCENT_LT=0.0" + "\n")
    new_file.write("PRIMER_INTERNAL_WT_HAIRPIN_TH=0.0" + "\n")
    new_file.write("PRIMER_INTERNAL_WT_LIBRARY_MISHYB=0.0" + "\n")
    new_file.write("PRIMER_INTERNAL_WT_NUM_NS=0.0" + "\n")
    new_file.write("PRIMER_INTERNAL_WT_SELF_ANY=0.0" + "\n")
    new_file.write("PRIMER_INTERNAL_WT_SELF_ANY_TH=0.0" + "\n")
    new_file.write("PRIMER_INTERNAL_WT_SELF_END=0.0" + "\n")
    new_file.write("PRIMER_INTERNAL_WT_SELF_END_TH=0.0" + "\n")
    new_file.write("PRIMER_INTERNAL_WT_SEQ_QUAL=0.0" + "\n")
    new_file.write("PRIMER_INTERNAL_WT_SIZE_GT=1.0" + "\n")
    new_file.write("PRIMER_INTERNAL_WT_SIZE_LT=1.0" + "\n")
    new_file.write("PRIMER_INTERNAL_WT_TM_GT=1.0" + "\n")
    new_file.write("PRIMER_INTERNAL_WT_TM_LT=1.0" + "\n")
    new_file.write("PRIMER_LIBERAL_BASE=1" + "\n")
    new_file.write("PRIMER_LIB_AMBIGUITY_CODES_CONSENSUS=0" + "\n")
    new_file.write("PRIMER_LOWERCASE_MASKING=0" + "\n")
    new_file.write("PRIMER_MAX_BOUND=110.0" + "\n")
    new_file.write("PRIMER_MAX_END_GC=5" + "\n")
    new_file.write("PRIMER_MAX_END_STABILITY=9.0" + "\n")
    new_file.write("PRIMER_MAX_GC=80.0" + "\n")
    new_file.write("PRIMER_MAX_HAIRPIN_TH=47.00" + "\n")
    new_file.write("PRIMER_MAX_LIBRARY_MISPRIMING=12.00" + "\n")
    new_file.write("PRIMER_MAX_NS_ACCEPTED=0" + "\n")
    new_file.write("PRIMER_MAX_POLY_X=5" + "\n")
    new_file.write("PRIMER_MAX_SELF_ANY=8.00" + "\n")
    new_file.write("PRIMER_MAX_SELF_ANY_TH=47.00" + "\n")
    new_file.write("PRIMER_MAX_SELF_END=3.00" + "\n")
    new_file.write("PRIMER_MAX_SELF_END_TH=47.00" + "\n")
    new_file.write("PRIMER_MAX_SIZE=27" + "\n")
    new_file.write("PRIMER_MAX_TEMPLATE_MISPRIMING=12.00" + "\n")
    new_file.write("PRIMER_MAX_TEMPLATE_MISPRIMING_TH=47.00" + "\n")
    new_file.write("PRIMER_MAX_TM=63.0" + "\n")
    new_file.write("PRIMER_MIN_3_PRIME_OVERLAP_OF_JUNCTION=4" + "\n")
    new_file.write("PRIMER_MIN_5_PRIME_OVERLAP_OF_JUNCTION=7" + "\n")
    new_file.write("PRIMER_MIN_BOUND=-10.0" + "\n")
    new_file.write("PRIMER_MIN_END_QUALITY=0" + "\n")
    new_file.write("PRIMER_MIN_GC=20.0" + "\n")
    new_file.write("PRIMER_MIN_LEFT_THREE_PRIME_DISTANCE=3" + "\n")
    new_file.write("PRIMER_MIN_QUALITY=0" + "\n")
    new_file.write("PRIMER_MIN_RIGHT_THREE_PRIME_DISTANCE=3" + "\n")
    new_file.write("PRIMER_MIN_SIZE=18" + "\n")
    new_file.write("PRIMER_MIN_TM=57.0" + "\n")
    new_file.write("PRIMER_NUM_RETURN=10" + "\n")
    new_file.write("PRIMER_OPT_BOUND=97.0" + "\n")
    new_file.write("PRIMER_OPT_GC_PERCENT=50.0" + "\n")
    new_file.write("PRIMER_OPT_SIZE=20" + "\n")
    new_file.write("PRIMER_OPT_TM=60.0" + "\n")
    new_file.write("PRIMER_OUTSIDE_PENALTY=0.0" + "\n")
    new_file.write("PRIMER_PAIR_MAX_COMPL_ANY=8.00" + "\n")
    new_file.write("PRIMER_PAIR_MAX_COMPL_ANY_TH=47.00" + "\n")
    new_file.write("PRIMER_PAIR_MAX_COMPL_END=3.00" + "\n")
    new_file.write("PRIMER_PAIR_MAX_COMPL_END_TH=47.00" + "\n")
    new_file.write("PRIMER_PAIR_MAX_DIFF_TM=100.0" + "\n")
    new_file.write("PRIMER_PAIR_MAX_LIBRARY_MISPRIMING=24.00" + "\n")
    new_file.write("PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING=24.00" + "\n")
    new_file.write("PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING_TH=47.00" + "\n")
    new_file.write("PRIMER_PAIR_WT_COMPL_ANY=0.0" + "\n")
    new_file.write("PRIMER_PAIR_WT_COMPL_ANY_TH=0.0" + "\n")
    new_file.write("PRIMER_PAIR_WT_COMPL_END=0.0" + "\n")
    new_file.write("PRIMER_PAIR_WT_COMPL_END_TH=0.0" + "\n")
    new_file.write("PRIMER_PAIR_WT_DIFF_TM=0.0" + "\n")
    new_file.write("PRIMER_PAIR_WT_IO_PENALTY=0.0" + "\n")
    new_file.write("PRIMER_PAIR_WT_LIBRARY_MISPRIMING=0.0" + "\n")
    new_file.write("PRIMER_PAIR_WT_PRODUCT_SIZE_GT=0.0" + "\n")
    new_file.write("PRIMER_PAIR_WT_PRODUCT_SIZE_LT=0.0" + "\n")
    new_file.write("PRIMER_PAIR_WT_PRODUCT_TM_GT=0.0" + "\n")
    new_file.write("PRIMER_PAIR_WT_PRODUCT_TM_LT=0.0" + "\n")
    new_file.write("PRIMER_PAIR_WT_PR_PENALTY=1.0" + "\n")
    new_file.write("PRIMER_PAIR_WT_TEMPLATE_MISPRIMING=0.0" + "\n")
    new_file.write("PRIMER_PAIR_WT_TEMPLATE_MISPRIMING_TH=0.0" + "\n")
    new_file.write("PRIMER_PICK_ANYWAY=0" + "\n")
    new_file.write("PRIMER_PICK_INTERNAL_OLIGO=0" + "\n")
    new_file.write("PRIMER_PICK_LEFT_PRIMER=1" + "\n")
    new_file.write("PRIMER_PICK_RIGHT_PRIMER=1" + "\n")
    new_file.write("PRIMER_PRODUCT_MAX_TM=1000000.0" + "\n")
    new_file.write("PRIMER_PRODUCT_MIN_TM=-1000000.0" + "\n")
    new_file.write("PRIMER_PRODUCT_OPT_SIZE=0" + "\n")
    new_file.write("PRIMER_PRODUCT_OPT_TM=0.0" + "\n")
    new_file.write("PRIMER_PRODUCT_SIZE_RANGE=600-800" + "\n")
    new_file.write("PRIMER_QUALITY_RANGE_MAX=100" + "\n")
    new_file.write("PRIMER_QUALITY_RANGE_MIN=0" + "\n")
    new_file.write("PRIMER_SALT_CORRECTIONS=1" + "\n")
    new_file.write("PRIMER_SALT_DIVALENT=1.5" + "\n")
    new_file.write("PRIMER_SALT_MONOVALENT=50.0" + "\n")
    new_file.write("PRIMER_SECONDARY_STRUCTURE_ALIGNMENT=1" + "\n")
    new_file.write("PRIMER_SEQUENCING_ACCURACY=20" + "\n")
    new_file.write("PRIMER_SEQUENCING_INTERVAL=250" + "\n")
    new_file.write("PRIMER_SEQUENCING_LEAD=50" + "\n")
    new_file.write("PRIMER_SEQUENCING_SPACING=500" + "\n")
    new_file.write("PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT=1" + "\n")
    new_file.write("PRIMER_THERMODYNAMIC_PARAMETERS_PATH=Add Path to the Primer3 /src/primer3_config/ folder" + "\n")
    new_file.write("PRIMER_THERMODYNAMIC_TEMPLATE_ALIGNMENT=0" + "\n")
    new_file.write("PRIMER_TM_FORMULA=1" + "\n")
    new_file.write("PRIMER_WT_BOUND_GT=0.0" + "\n")
    new_file.write("PRIMER_WT_BOUND_LT=0.0" + "\n")
    new_file.write("PRIMER_WT_END_QUAL=0.0" + "\n")
    new_file.write("PRIMER_WT_END_STABILITY=0.0" + "\n")
    new_file.write("PRIMER_WT_GC_PERCENT_GT=0.0" + "\n")
    new_file.write("PRIMER_WT_GC_PERCENT_LT=0.0" + "\n")
    new_file.write("PRIMER_WT_HAIRPIN_TH=0.0" + "\n")
    new_file.write("PRIMER_WT_LIBRARY_MISPRIMING=0.0" + "\n")
    new_file.write("PRIMER_WT_NUM_NS=0.0" + "\n")
    new_file.write("PRIMER_WT_POS_PENALTY=1.0" + "\n")
    new_file.write("PRIMER_WT_SELF_ANY=0.0" + "\n")
    new_file.write("PRIMER_WT_SELF_ANY_TH=0.0" + "\n")
    new_file.write("PRIMER_WT_SELF_END=0.0" + "\n")
    new_file.write("PRIMER_WT_SELF_END_TH=0.0" + "\n")
    new_file.write("PRIMER_WT_SEQ_QUAL=0.0" + "\n")
    new_file.write("PRIMER_WT_SIZE_GT=1.0" + "\n")
    new_file.write("PRIMER_WT_SIZE_LT=1.0" + "\n")
    new_file.write("PRIMER_WT_TEMPLATE_MISPRIMING=0.0" + "\n")
    new_file.write("PRIMER_WT_TEMPLATE_MISPRIMING_TH=0.0" + "\n")
    new_file.write("PRIMER_WT_TM_GT=1.0" + "\n")
    new_file.write("PRIMER_WT_TM_LT=1.0" + "\n")
    new_file.write("SEQUENCE_EXCLUDED_REGION=400,200" + "\n")
    new_file.write("SEQUENCE_FORCE_LEFT_END=-1000000" + "\n")
    new_file.write("SEQUENCE_FORCE_LEFT_START=-1000000" + "\n")
    new_file.write("SEQUENCE_FORCE_RIGHT_END=-1000000" + "\n")
    new_file.write("SEQUENCE_FORCE_RIGHT_START=-1000000" + "\n")
    new_file.write("SEQUENCE_START_CODON_POSITION=-2000000" + "\n")
    new_file.write("SEQUENCE_START_CODON_SEQUENCE=ATG" + "\n")
    new_file.write("=\n")
    new_file.close()
"""   
Execute the primer3_core command on the input file you created.
"""

primer3_command = 'primer3_core < primer.txt'
output = subprocess.check_output(primer3_command, shell=True, encoding = 'utf-8')

""" 
Read in and parse the output to identify the sequence of best two primers. Print these out to the user
"""

primer_output = output.split('\n')
best_primers = []
for i, line in enumerate(primer_output):
    if line.startswith("PRIMER_LEFT_0_SEQUENCE="):
        best_primers.append(line.split('=')[1])
    elif line.startswith("PRIMER_RIGHT_0_SEQUENCE="):
        best_primers.append(line.split('=')[1])
    if len(best_primers) == 2:
        break

for i, primer in enumerate(best_primers):
    print(f"{primer}")

        