#!/usr/bin/env python3
"""
rhodopsin_neighborhood.py

This script analyzes DRAM's output CSV table to identify rhodopsin genes and their genomic neighborhoods,
detecting features like multiple rhodopsins, carotenoid gene clusters, and flotillin associations.
"""

import pandas as pd
import argparse
import sys
from typing import Dict, List, Any
from concurrent.futures import ThreadPoolExecutor, as_completed

def find_rhodopsin_genes(df: pd.DataFrame) -> pd.DataFrame:
    """
    Find rows where pfam_hits contain the string 'rhodopsin'
    """
    rhodopsin_mask = df['pfam_hits'].str.contains('rhodopsin', case=False, na=False)
    return df[rhodopsin_mask].copy()

def count_rhodopsin_per_sag(rhodopsin_df: pd.DataFrame) -> Dict[str, int]:
    """
    Count the number of rhodopsin genes in each SAG (fasta)
    Returns dictionary with fasta as key and count as value
    """
    fasta_counts = rhodopsin_df['fasta'].value_counts()
    return fasta_counts.to_dict()

def get_flanking_genes(df: pd.DataFrame, rhodopsin_row: pd.Series, window_size: int = 5) -> Dict[str, List]:
    """
    Extract the 5 flanking genes from each side of the rhodopsin gene, including the rhodopsin gene itself
    Returns dictionary with annot_multisrc as keys and [annot_multisrc, hmmsearch_gene, start_position, end_position, strandedness] as values
    """
    scaffold = rhodopsin_row['scaffold']
    gene_position = rhodopsin_row['gene_position']
    
    # Get all genes on the same scaffold
    scaffold_genes = df[df['scaffold'] == scaffold].sort_values('gene_position')
    
    # Find the index of the rhodopsin gene
    rhodopsin_idx = scaffold_genes[scaffold_genes['...1'] == rhodopsin_row['...1']].index
    if len(rhodopsin_idx) == 0:
        return {}
    
    rhodopsin_idx = rhodopsin_idx[0]
    scaffold_indices = scaffold_genes.index.tolist()
    rhodopsin_pos = scaffold_indices.index(rhodopsin_idx)
    
    # Calculate flanking gene indices
    start_idx = max(0, rhodopsin_pos - window_size)
    end_idx = min(len(scaffold_indices), rhodopsin_pos + window_size + 1)
    
    flanking_genes = {}
    for i in range(start_idx, end_idx):
        gene_row = scaffold_genes.loc[scaffold_indices[i]]
        
        # Get strandedness from the 'strandedness' column
        strandedness = int(gene_row.get('strandedness', 0))  # Default to 0 if not found
        
        # Determine the gene position relative to the rhodopsin gene
        if i < rhodopsin_pos:
            pos = i - rhodopsin_pos  # Negative for upstream genes
            annot_key = f"annot_multisrc of gene {pos}"
        elif i > rhodopsin_pos:
            pos = i - rhodopsin_pos  # Positive for downstream genes
            annot_key = f"annot_multisrc of gene +{pos}"
        else:
            annot_key = "annot_multisrc of rhodopsin gene"
            
        flanking_genes[annot_key] = [
            str(gene_row['annot_multisrc']) if pd.notna(gene_row['annot_multisrc']) else "None",
            str(gene_row['hmmsearch_gene']) if pd.notna(gene_row['hmmsearch_gene']) else "None",
            int(gene_row['start_position']),
            int(gene_row['end_position']),
            int(strandedness)
        ]
    
    return flanking_genes

def count_carotenoid_genes(flanking_genes: Dict[str, List]) -> int:
    """
    Count detection times of carotenoid-related strings in flanking genes
    Uses exact case-sensitive matching
    """
    carotenoid_terms = ["CrtB", "CrtE", "CrtI", "CrtY", "brp/blh"]
    
    count = 0
    for gene_info in flanking_genes.values():
        hmmsearch_gene = gene_info[1] if len(gene_info) > 0 else ""
        
        # Check for exact case-sensitive matches
        gene_value = str(hmmsearch_gene).strip()
        for term in carotenoid_terms:
            if gene_value == term:  # Exact case-sensitive match
                count += 1
                break  # Count each gene only once
    
    return count

def detect_flotillin_adjacent(df: pd.DataFrame, rhodopsin_row: pd.Series) -> bool:
    """
    Check if flotillin is detected in the most adjacent gene (upstream or downstream)
    """
    scaffold = rhodopsin_row['scaffold']
    gene_position = rhodopsin_row['gene_position']
    
    # Get all genes on the same scaffold
    scaffold_genes = df[df['scaffold'] == scaffold].sort_values('gene_position')
    
    # Find the index of the rhodopsin gene
    rhodopsin_idx = scaffold_genes[scaffold_genes['...1'] == rhodopsin_row['...1']].index
    if len(rhodopsin_idx) == 0:
        return False
    
    rhodopsin_idx = rhodopsin_idx[0]
    scaffold_indices = scaffold_genes.index.tolist()
    rhodopsin_pos = scaffold_indices.index(rhodopsin_idx)
    
    # Check immediate upstream neighbor
    if rhodopsin_pos > 0:
        upstream_gene = scaffold_genes.loc[scaffold_indices[rhodopsin_pos - 1]]
        if pd.notna(upstream_gene['hmmsearch_gene']) and 'flotillin' in str(upstream_gene['hmmsearch_gene']).lower():
            return True
    
    # Check immediate downstream neighbor
    if rhodopsin_pos < len(scaffold_indices) - 1:
        downstream_gene = scaffold_genes.loc[scaffold_indices[rhodopsin_pos + 1]]
        if pd.notna(downstream_gene['hmmsearch_gene']) and 'flotillin' in str(downstream_gene['hmmsearch_gene']).lower():
            return True
    
    return False

def process_rhodopsin_gene(rhodopsin_row: pd.Series, df: pd.DataFrame, rhodopsin_counts_per_sag: Dict[str, int]) -> Dict[str, Any]:
    """
    Process a single rhodopsin gene and return its analysis results.
    This function is designed to be run in parallel.
    """
    print(f"Processing rhodopsin gene: {rhodopsin_row['...1']}")
    
    # (3) Get flanking genes
    flanking_genes = get_flanking_genes(df, rhodopsin_row)
    
    # (4) Count carotenoid genes in flanking regions
    carotenoid_count = count_carotenoid_genes(flanking_genes)
    carotenoid_cluster = 1 if carotenoid_count >= 2 else 0
    
    # (5) Check for flotillin in adjacent genes
    flotillin_detected = detect_flotillin_adjacent(df, rhodopsin_row)
    
    # Check if there are 5 flanking genes on both sides (10 total) and neither carotenoid nor flotillin
    neither = 0
    if len(flanking_genes) == 10:  # 5 genes on each side
        if not carotenoid_cluster and not flotillin_detected:
            neither = 1
    
    # Prepare and return output row
    return {
        'rhodopsin_gene': rhodopsin_row['...1'],
        'rhodopsin_scaffold': rhodopsin_row['scaffold'],
        'sag': rhodopsin_row['fasta'],
        'rhodopsin_count_in_sag': rhodopsin_counts_per_sag.get(rhodopsin_row['fasta'], 0),
        'carotenoid_gene_cluster': carotenoid_cluster,
        'flotillin': 1 if flotillin_detected else 0,
        'neither': neither,
        'flanking_genes': flanking_genes
    }

def main():
    parser = argparse.ArgumentParser(description='Analyze rhodopsin genes and their genomic neighborhoods')
    parser.add_argument('input_file', help='Input CSV file containing gene annotations')
    parser.add_argument('-o', '--output', default='rhodopsin_analysis.csv', help='Output CSV file name')
    parser.add_argument('--max_workers', type=int, default=None, 
                       help='Maximum number of worker threads to use')
    
    args = parser.parse_args()
    
    try:
        # Read input CSV
        print(f"Reading input file: {args.input_file}")
        df = pd.read_csv(args.input_file)
        
        # Check required columns
        required_columns = ['...1', 'fasta', 'scaffold', 'gene_position', 'start_position', 'end_position', 'kegg_hit', 'pfam_hits', 'annot_multisrc']
        missing_columns = [col for col in required_columns if col not in df.columns]
        if missing_columns:
            print(f"Error: Missing required columns: {missing_columns}")
            sys.exit(1)
        
        # (1) Find rhodopsin genes
        print("Finding rhodopsin genes...")
        rhodopsin_df = find_rhodopsin_genes(df)
        
        if rhodopsin_df.empty:
            print("No rhodopsin genes found in the dataset.")
            return
        
        print(f"Found {len(rhodopsin_df)} rhodopsin genes")
        
        # (2) Count rhodopsin genes per SAG
        print("Counting rhodopsin genes per SAG...")
        rhodopsin_counts_per_sag = count_rhodopsin_per_sag(rhodopsin_df)

        counts_series = pd.Series(rhodopsin_counts_per_sag)
        
        # Process rhodopsin genes in parallel
        print(f"Processing {len(rhodopsin_df)} rhodopsin genes with {args.max_workers or 'default'} workers...")
        output_data = []
        
        # Use ThreadPoolExecutor for parallel processing
        with ThreadPoolExecutor(max_workers=args.max_workers) as executor:
            # Submit all tasks
            future_to_row = {
                executor.submit(process_rhodopsin_gene, row, df, rhodopsin_counts_per_sag): idx 
                for idx, row in rhodopsin_df.iterrows()
            }
            
            # Process results as they complete
            for future in as_completed(future_to_row):
                try:
                    result = future.result()
                    output_data.append(result)
                except Exception as exc:
                    idx = future_to_row[future]
                    print(f"Error processing gene at index {idx}: {exc}")
        
        # Create output DataFrame
        output_df = pd.DataFrame(output_data)
        
        # Save results
        print(f"Saving results to: {args.output}")
        output_df.to_csv(args.output, index=False)
        
        # Print summary
        print("\n=== Analysis Summary ===")
        print(f"Total rhodopsin genes analyzed: {len(output_df)}")
        print(f"SAGs with multiple rhodopsins: {(counts_series > 1).sum()}")
        print(f"SAGs with single rhodopsins: {(counts_series == 1).sum()}")
        print(f"Genes with carotenoid gene clusters: {output_df['carotenoid_gene_cluster'].sum()}")
        print(f"Genes with adjacent flotillin: {output_df['flotillin'].sum()}")
        print(f"Rhodopsin count distribution:")
        print(f"  Min: {counts_series.min()}")
        print(f"  Max: {counts_series.max()}")
        
    except FileNotFoundError:
        print(f"Error: Input file '{args.input_file}' not found.")
        sys.exit(1)
    except Exception as e:
        print(f"Error during analysis: {str(e)}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()
