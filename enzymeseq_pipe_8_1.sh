#!/bin/bash
set -euo pipefail
topdir=$(pwd)
account="b1042"
queue="genomics"
DEBUG="$topdir/debug"
total_dir="$topdir/total_output"
polyA_dir="$topdir/polyA_output"
micropolyA_dir="$topdir/micropolyA_output"
OUTDIR="$topdir/micropolyA_dedup"
mkdir -p "$OUTDIR" "$DEBUG"

for sample_dir in "$micropolyA_dir"/*; do
    [[ -d "$sample_dir" ]] || continue
    sample=$(basename "$sample_dir")
    echo "Submitting job for sample: $sample"
    mkdir -p "$OUTDIR/$sample"
    export sample
	
 sbatch --export=ALL,MICROPOLYA_DIR="$micropolyA_dir",POLYA_DIR="$polyA_dir",OUTDIR="$OUTDIR",DEBUG="$DEBUG" <<EOF
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

export CONDA_BACKUP_ADDR2LINE=""
set +u
source ~/.bashrc
conda activate seqkit
set -u

for micropolyA_file in "\$MICROPOLYA_DIR/\$sample/"*_utr.fasta; do
    [[ -s "\$micropolyA_file" ]] || continue

    gene=\$(basename "\$micropolyA_file" _utr.fasta)
    outdir="\$OUTDIR/\$sample"
    mkdir -p "\$outdir"
    dedup_file="\$outdir/\${gene}_micropolyA_dedup.fasta"

    # Check if a polyA file exists
    polyA_file="\$POLYA_DIR/\$sample/polyA_matched_reads/\${gene}_utr.fasta"

    if [[ -s "\$polyA_file" ]]; then
        echo "Deduplicating \$gene using polyA sequences"
        seqkit seq -n "\$polyA_file" > "\$outdir/\${gene}_polyA_ids.txt"
        seqkit grep -v -r -n -f "\$outdir/\${gene}_polyA_ids.txt" "\$micropolyA_file" > "\$dedup_file"
    else
        echo "No polyA file for \$gene, keeping all micropolyA sequences"
        cp "\$micropolyA_file" "\$dedup_file"
    fi

    [[ -s "\$dedup_file" ]] || echo "WARNING: deduplicated file is empty for \$gene"
done
EOF

done