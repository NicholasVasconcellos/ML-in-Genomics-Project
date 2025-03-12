#!/usr/bin/env python

import sys
import random
from Bio import SeqIO

def random_deletion(input_file, output_file, num_to_delete):
    try:
        to_delete = int(num_to_delete)
    except ValueError:
        sys.exit(f"Error: Invalid number of sequences to delete: {num_to_delete}")
    
    # Read all sequences
    records = list(SeqIO.parse(input_file, 'fasta'))
    if len(records) <= 1:
        sys.exit('Error: MSA contains only query sequence')
    
    # Number of sequences to delete (ensure we keep at least 1 non-query sequence)
    to_delete = min(to_delete, len(records) - 2)
    
    if to_delete <= 0:
        # If parameter is invalid, just copy the original file
        with open(input_file, 'r') as f_in, open(output_file, 'w') as f_out:
            f_out.write(f_in.read())
    else:
        # Always keep query sequence (first record)
        to_keep = [0]
        # Randomly select sequences to keep
        to_keep.extend(random.sample(range(1, len(records)), len(records) - 1 - to_delete))
        to_keep.sort()
        
        # Write kept sequences to new file
        with open(output_file, 'w') as f_out:
            for idx in to_keep:
                f_out.write(f'>{records[idx].id}\n{records[idx].seq}\n')

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print(f"Usage: {sys.argv[0]} input_file output_file num_to_delete")
        sys.exit(1)
    
    random_deletion(sys.argv[1], sys.argv[2], sys.argv[3])