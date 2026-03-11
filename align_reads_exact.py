#!/usr/bin/env python3

from Bio import SeqIO
import argparse
import numpy as np
from multiprocessing import Pool
import functools
import os
import re

# ----------------------------
# Utilities
# ----------------------------

def calculate_mismatches(seq1, seq2):
    """Count mismatches using numpy"""
    if len(seq1) == 0:
        return 0
    arr1 = np.frombuffer(seq1.encode(), dtype=np.uint8)
    arr2 = np.frombuffer(seq2.encode(), dtype=np.uint8)
    return np.sum(arr1 != arr2)

def align_read_full(read_seq, ref_seq, max_mismatch_rate=0.0):
    """Check if full read matches reference 3' end"""
    rlen = len(read_seq)
    if rlen == 0:
        return False
    ref_len = len(ref_seq)
    if rlen > ref_len:
        return False
    ref_segment = ref_seq[-rlen:]
    mismatches = calculate_mismatches(read_seq, ref_segment)
    return (mismatches / rlen) <= max_mismatch_rate

def process_batch(batch, refs, gene_id, max_mismatch_rate):
    """Align a batch of reads"""
    results = []
    for record_id, read_seq in batch:
        if gene_id not in refs:
            results.append((record_id, read_seq, False))
        else:
            aligned = align_read_full(read_seq, refs[gene_id], max_mismatch_rate)
            results.append((record_id, read_seq, aligned))
    return results



def main(reads_file, refs_file, aligned_fasta, unaligned_fasta,
         gene_id, max_mismatch_rate=0.0, n_processes=4, batch_size=1000, count_empty_as_one=False):
    
   
    refs = {}
    for rec in SeqIO.parse(refs_file, "fasta"):
        gene_id_match = re.match(r'(VACWR\d+)', rec.id)
        if gene_id_match:
            refs[gene_id_match.group(1)] = str(rec.seq).upper()
    print(f"Loaded {len(refs)} reference sequences")

    if gene_id not in refs:
        print(f"WARNING: Gene {gene_id} not found in reference")

    reads_data = [(r.id, str(r.seq).upper()) for r in SeqIO.parse(reads_file, "fasta")]
    print(f"Loaded {len(reads_data)} reads from {reads_file} assigned to gene {gene_id}")

    batches = [reads_data[i:i+batch_size] for i in range(0, len(reads_data), batch_size)]
    process_func = functools.partial(process_batch, refs=refs, gene_id=gene_id, max_mismatch_rate=max_mismatch_rate)

    aligned_count = 0
    unaligned_count = 0
    empty_count = 0

    os.makedirs(os.path.dirname(aligned_fasta), exist_ok=True)

    with open(aligned_fasta, 'w') as af, open(unaligned_fasta, 'w') as uf:
        with Pool(n_processes) as pool:
            for batch_results in pool.imap(process_func, batches):
                for record_id, read_seq, aligned in batch_results:
                    # Count empty reads
                    if len(read_seq) == 0:
                        empty_count += 1
                        if count_empty_as_one:
                            aligned = True
                            read_seq = "N"
                    if aligned:
                        af.write(f">{record_id}\n{read_seq}\n")
                        aligned_count += 1
                    else:
                        uf.write(f">{record_id}\n{read_seq}\n")
                        unaligned_count += 1

    print(f"Alignment finished:")
    print(f"  Aligned reads: {aligned_count}")
    print(f"  Unaligned reads: {unaligned_count}")
    print(f"  Empty reads: {empty_count}")
    print(f"Aligned FASTA: {aligned_fasta}")
    print(f"Unaligned FASTA: {unaligned_fasta}")



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Full-read 3'-anchored alignment for VACV")
    parser.add_argument("reads", help="Input FASTA reads")
    parser.add_argument("refs", help="Reference FASTA")
    parser.add_argument("aligned_fasta", help="Output aligned FASTA (_highscore_aligned_seqs.fa)")
    parser.add_argument("unaligned_fasta", help="Output unaligned FASTA (_unaligned.fa)")
    parser.add_argument("--gene", required=True, help="Gene ID for this file (from filename)")
    parser.add_argument("--mismatch-rate", type=float, default=0.0, help="Max mismatch rate (default 0.0)")
    parser.add_argument("--processes", type=int, default=4, help="Number of parallel processes")
    parser.add_argument("--batch-size", type=int, default=1000, help="Batch size for multiprocessing")
    parser.add_argument("--count-empty-as-one", action="store_true", help="Optionally count empty reads as aligned placeholder")
    args = parser.parse_args()

    main(args.reads, args.refs, args.aligned_fasta, args.unaligned_fasta,
         args.gene, args.mismatch_rate, args.processes, args.batch_size, args.count_empty_as_one)
