#!/usr/bin/env python

import sys

def taxonomic_deletion(input_file, output_file, remove_file):
    # Read list of sequences to remove
    with open(remove_file, 'r') as f:
        remove = set(line.strip() for line in f)
    
    # Process the MSA file
    keep_next = True
    with open(input_file, 'r') as f_in, open(output_file, 'w') as f_out:
        for line in f_in:
            if line.startswith('>'):
                header = line.strip().split()[0][1:]
                keep_next = header not in remove
                if keep_next:
                    f_out.write(line)
            elif keep_next:
                f_out.write(line)

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print(f"Usage: {sys.argv[0]} input_file output_file remove_file")
        sys.exit(1)
    
    taxonomic_deletion(sys.argv[1], sys.argv[2], sys.argv[3])