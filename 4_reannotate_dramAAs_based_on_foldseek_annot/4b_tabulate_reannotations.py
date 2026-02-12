#!/usr/bin/env python3
"""
Script to extract information from proteins with 'Rank E' (or 'rank: G') in their description
and create a CSV table with gene_id, description, descr_score, best_score, and pfam_acc.

This version uses Biopython for FASTA parsing.

USAGE:

python 4b_tabulate_reannotations.py reannot_Erank_proteins_from_rhodopsin_contigs.fasta -o rankE_reannotations.csv -r G

"""

import re
import csv
from Bio import SeqIO


def parse_header(header):
    """
    Parse the FASTA header to extract relevant information.
    
    Expected format:
    >AG-538-A02_AG-538-A02_contigs_AG-538-A02_NODE_4_3 rank: G; Reverse transcriptase domain-containing protein bit-148/440 
    
    Returns:
        dict with keys: gene_id, rank, description, descr_score, best_score, pfam_acc
    """
    # Initialize return dictionary
    result = {
        'gene_id': '',
        'rank': '',
        'description': '',
        'descr_score': '',
        'best_score': '',
        'pfam_acc': ''
    }
    
    # Split header into parts
    parts = header.split(None, 1)  # Split on first whitespace
    if len(parts) < 2:
        return None
    
    # Extract gene_id (remove the prefix portion to match CSV format)
    full_id = parts[0]
    # Extract the simplified gene_id format (e.g., AG-538-A02_NODE_4_3)
    match = re.search(r'([A-Z]+-\d+-[A-Z]\d+)_.+_([A-Z]+-\d+-[A-Z]\d+_NODE_\d+_\d+)', full_id)
    if match:
        result['gene_id'] = f"{match.group(1)}_{match.group(2).split('_', 1)[1]}"
    else:
        # Alternative pattern for SCGC entries
        match = re.search(r'(SCGC_[A-Z]+-\d+-[A-Z]\d+)_(.+)', full_id)
        if match:
            result['gene_id'] = f"{match.group(1)}_{match.group(2)}"
        else:
            result['gene_id'] = full_id
    
    # Parse the rest of the header
    rest = parts[1]
    
    # Extract rank
    rank_match = re.search(r'rank:\s*([A-Z])', rest)
    if rank_match:
        result['rank'] = rank_match.group(1)
    
    # Extract description and bit scores
    # Pattern: description bit-descr_score/best_score OR just description
    desc_pattern = r'rank:\s*[A-Z];\s*(.+?)(?:\s+bit-(\d+)/(\d+))?(?:\s+\[([^\]]+)\])?(?:\s+\(db=.*\))?\s*$'
    desc_match = re.search(desc_pattern, rest)
    
    if desc_match:
        result['description'] = desc_match.group(1).strip()
        if desc_match.group(2):
            result['descr_score'] = desc_match.group(2)
        if desc_match.group(3):
            result['best_score'] = desc_match.group(3)
        if desc_match.group(4):
            result['pfam_acc'] = desc_match.group(4)
    
    return result


def extract_rank_e_proteins(fasta_file, output_csv, target_rank='G'):
    """
    Extract proteins with specified rank from FASTA file and write to CSV.
    
    Args:
        fasta_file: Path to input FASTA file
        output_csv: Path to output CSV file
        target_rank: The rank letter to filter for (default 'G', but can be 'E')
    """
    results = []
    
    # Parse FASTA file using Biopython
    with open(fasta_file, 'r') as f:
        for record in SeqIO.parse(f, 'fasta'):
            # Parse header
            parsed = parse_header(record.description)
            
            if parsed and parsed['rank'] == target_rank:
                results.append({
                    'gene_id': parsed['gene_id'],
                    'description': parsed['description'],
                    'descr_score': parsed['descr_score'],
                    'best_score': parsed['best_score'],
                    'pfam_acc': parsed['pfam_acc']
                })
    
    # Keep original FASTA order (no sorting)
    
    # Write to CSV
    with open(output_csv, 'w', newline='') as csvfile:
        fieldnames = ['gene_id', 'description', 'descr_score', 'best_score', 'pfam_acc']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        
        writer.writeheader()
        for row in results:
            writer.writerow(row)
    
    print(f"Extracted {len(results)} proteins with rank '{target_rank}'")
    print(f"Output written to: {output_csv}")


if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(
        description='Extract Rank E (or other rank) proteins from FASTA and create CSV table (Biopython version)'
    )
    parser.add_argument(
        'input_fasta',
        help='Input FASTA file'
    )
    parser.add_argument(
        '-o', '--output',
        default='rank_e_proteins.csv',
        help='Output CSV file (default: rank_e_proteins.csv)'
    )
    parser.add_argument(
        '-r', '--rank',
        default='G',
        help='Rank letter to filter for (default: G). Note: The file is labeled "Erank" but contains "G" ranks'
    )
    
    args = parser.parse_args()
    
    extract_rank_e_proteins(args.input_fasta, args.output, args.rank)
