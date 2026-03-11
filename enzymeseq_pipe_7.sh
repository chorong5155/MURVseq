#!/bin/bash
set -euo pipefail
topdir=$(pwd)
account="b1042"
queue="genomics"
DEBUG="$topdir/debug"
total_dir="$topdir/total_output"
aligned_dir="$topdir/aligned_UTR"
OUTDIR="$topdir/un_natural_output"
mkdir -p "$OUTDIR" "$DEBUG"
export topdir account queue DEBUG total_dir aligned_dir OUTDIR


for sample_dir in "$aligned_dir"/*; do
    [[ -d "$sample_dir" ]] || continue
    sample=$(basename "$sample_dir")
    echo "Submitting job for sample: $sample"

   
    mkdir -p "$OUTDIR/$sample"
    export sample

    sbatch --export=ALL <<EOF
#!/bin/bash -l
#SBATCH -A $account
#SBATCH -p $queue
#SBATCH -o $DEBUG/seqkit_${sample}.%j.out
#SBATCH -e $DEBUG/seqkit_${sample}.%j.err
#SBATCH -t 10:00:00
#SBATCH -c 4
#SBATCH --mem-per-cpu=10G
#SBATCH -J seqkit_${sample}

set -euo pipefail

# Backup for conda
export CONDA_BACKUP_ADDR2LINE=""
set +u
source ~/.bashrc
conda activate seqkit
set -u

# Directories
highscore_dir="\$aligned_dir/\$sample"
full_utr_dir="\$total_dir/\$sample"
sample_outdir="\$OUTDIR/\$sample"
mkdir -p "\$sample_outdir/totalUTRs_cleaned" \
         "\$sample_outdir/naturalUTRs" \
         "\$sample_outdir/unnaturalUTRs"

shopt -s nullglob
highscore_files=(\$highscore_dir/*_utr_highscore_aligned_seqs.fa)


for highscore_fa in "\${highscore_files[@]}"; do
    [[ -f "\$highscore_fa" ]] || continue
    gene_name=\$(basename "\$highscore_fa" "_utr_highscore_aligned_seqs.fa")
    full_utr_fa="\$full_utr_dir/\${gene_name}_utr.fasta"

    # Output files
    natural_id_file_clean="\$sample_outdir/\${gene_name}_naturalUTRs_ids.txt.clean"
    full_utr_fa_clean="\$sample_outdir/totalUTRs_cleaned/\${gene_name}_total_utr_cleaned.fasta"
    natural_out="\$sample_outdir/naturalUTRs/\${gene_name}_naturalUTRs.fasta"
    unnatural_out="\$sample_outdir/unnaturalUTRs/\${gene_name}_unnaturalUTRs.fasta"

    # Skip if highscore file is empty
    if [[ ! -s "\$highscore_fa" ]]; then
        echo "Highscore file \$highscore_fa is empty, creating empty output files and skipping."
        touch "\$natural_out" "\$unnatural_out"
        continue
    fi

    # Skip if full UTR file does not exist
    [[ -f "\$full_utr_fa" ]] || continue

    # Extract IDs from highscore FASTA
    # Extract only the part before the first space
seqkit seq -n "\$highscore_fa" | sort -u > "\$natural_id_file_clean"
    num_ids=\$(wc -l < "\$natural_id_file_clean")
if [[ "\$num_ids" -eq 0 ]]; then
    echo "WARNING: No IDs extracted from \$highscore_fa"
    touch "\$natural_out" "\$unnatural_out"
    continue
fi
seqkit replace -p ' .*' -r '' "\$full_utr_fa" -o "\$full_utr_fa_clean"
num_full=\$(grep -c '^>' "\$full_utr_fa_clean")
if [[ "\$num_full" -eq 0 ]]; then
    echo "WARNING: Full UTR file \$full_utr_fa_clean has no sequences after cleaning headers"
    touch "\$natural_out" "\$unnatural_out"
    continue
fi
    # Extract natural UTRs
    seqkit grep -n -f "\$natural_id_file_clean" "\$full_utr_fa_clean" > "\$natural_out"

    # Extract unnatural UTRs
    seqkit grep -n -v -f "\$natural_id_file_clean" "\$full_utr_fa_clean" > "\$unnatural_out"
num_unnatural=\$(grep -c '^>' "\$unnatural_out")
if [[ "\$num_unnatural" -eq 0 ]]; then
    echo "WARNING: No unnatural UTRs found for gene \$gene_name in sample \$sample"
fi
    # Cleanup
    rm -f "\$natural_id_file_clean"
done

echo "Sample \$sample processing complete."
EOF

done

