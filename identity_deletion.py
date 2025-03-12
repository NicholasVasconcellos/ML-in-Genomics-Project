#!/usr/bin/env python

import sys
from Bio import SeqIO

def calc_identity(seq1, seq2):
    """Calculate sequence identity between two sequences"""
    # Remove gaps for identity calculation
    seq1_nogap = seq1.replace('-', '')
    seq2_aligned = ''
    idx = 0
    for char in seq2:
        if char != '-':
            if idx < len(seq1_nogap):
                seq2_aligned += char
            idx += 1
    
    # Calculate identity over aligned region
    matches = sum(1 for a, b in zip(seq1_nogap[:len(seq2_aligned)], seq2_aligned) if a == b)
    return matches / len(seq1_nogap[:len(seq2_aligned)]) if len(seq1_nogap[:len(seq2_aligned)]) > 0 else 0

def identity_deletion(input_file, output_file, threshold):
    threshold = float(threshold)
    
    # Read all sequences
    records = list(SeqIO.parse(input_file, 'fasta'))
    if not records:
        sys.exit('Error: No sequences found in MSA')
    
    # Extract query sequence
    query_seq = str(records[0].seq).replace('-', '')
    
    # Write filtered MSA
    with open(output_file, 'w') as f_out:
        # Always keep the query sequence
        f_out.write(f'>{records[0].id}\n{records[0].seq}\n')
        
        # Filter other sequences based on identity
        for record in records[1:]:
            seq_identity = calc_identity(query_seq, str(record.seq))
            # Keep sequences with identity LESS than threshold
            if seq_identity < threshold:
                f_out.write(f'>{record.id}\n{record.seq}\n')

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print(f"Usage: {sys.argv[0]} input_file output_file threshold")
        sys.exit(1)
    
    identity_deletion(sys.argv[1], sys.argv[2], sys.argv[3])