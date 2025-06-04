from Bio import SeqIO
#import pandas as pd
import sys
import os
import re

#remove sequences with more than 10% N or ?
def filter_and_convert_fasta(fasta_file, output_fa=None):
    with open(output_fa, "w") as output_handle:
        for record in SeqIO.parse(fasta_file, "fasta"):
            sequence = str(record.seq)
            invalid_count = sequence.upper().count('N') + sequence.upper().count('?') + sequence.upper().count('-')
            if invalid_count / len(sequence) <= 0.1:  # Keep sequences with <= 10% N, ?, or -
                output_handle.write(f">{record.id}\n{sequence}\n")

# Usage
if __name__ == "__main__":
    fasta_file = sys.argv[1]  # Replace with your FASTA file
    fname = os.path.basename(fasta_file)
    fname = re.sub(r'\.fa(sta)?$', '', fname)  # Handle .fa or .fasta extensions
    output_fa = f'{fname}_del.fa'  # CSV output file

    filter_and_convert_fasta(fasta_file, output_fa)
