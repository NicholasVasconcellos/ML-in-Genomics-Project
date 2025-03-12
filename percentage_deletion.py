#!/usr/bin/env python

import sys
import random
from Bio import SeqIO

def percentage_deletion(input_file, output_file, percentage):
    try:
        keep_percent = float(percentage) / 100.0
        if keep_percent < 0 or keep_percent > 1:
            raise ValueError("Percentage must be between 0 and 100")
    except ValueError as e:
        sys.exit(f"Error: Invalid percentage value: {e}")
    
    # Read all sequences
    records = list(SeqIO.parse(input_file, 'fasta'))
    if len(records) <= 1:
        sys.exit('Error: MSA contains only query sequence')
    
    # Calculate number to keep based on percentage
    to_keep_count = max(1, int((len(records) - 1) * keep_percent)) + 1  # +1 for query
    
    # Always keep query sequence (first record)
    to_keep = [0]
    # Randomly select sequences to keep
    to_keep.extend(random.sample(range(1, len(records)), min(to_keep_count - 1, len(records) - 1)))
    to_keep.sort()
    
    # Write kept sequences to new file
    with open(output_file, 'w') as f_out:
        for idx in to_keep:
            f_out.write(f'>{records[idx].id}\n{records[idx].seq}\n')

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print(f"Usage: {sys.argv[0]} input_file output_file percentage")
        sys.exit(1)
    
    percentage_deletion(sys.argv[1], sys.argv[2], sys.argv[3])