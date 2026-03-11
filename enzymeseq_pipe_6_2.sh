#!/bin/bash

topdir=$(pwd)
account="b1042"
queue="genomicslong"
DEBUG="${topdir}/debug"
refDir="/projects/b1154/referenceGenomes"
genomeID="hg38+VacV"
GENE_REF_FA="${refDir}/${genomeID}/VacV_150utr_flanks.fa"
total_dir="${topdir}/total_output"
dedup_dir="${topdir}/dedup_output"
aligned_output="${topdir}/aligned_UTR_2"
PY_SCRIPT="/projects/b1154/cpp4747/scripts/align_reads_exact.py"

mkdir -p "$DEBUG" "$dedup_dir" "$aligned_output"

for sample_dir in "$total_dir"/*; do
    [[ -d "$sample_dir" ]] || continue
    sample=$(basename "$sample_dir")

    echo "Submitting job for sample: $sample"

    mkdir -p "$dedup_dir/$sample"
    mkdir -p "$aligned_output/$sample"

    sbatch  --export=ALL,sample_dir="$sample_dir",sample="$sample",dedup_dir="$dedup_dir",aligned_output="$aligned_output",GENE_REF_FA="$GENE_REF_FA",PY_SCRIPT="$PY_SCRIPT" <<EOF
#!/bin/bash -l
#SBATCH -A $account
#SBATCH --partition=$queue
#SBATCH -J aligned_UTR_${sample}
#SBATCH -o $DEBUG/aligned_UTR_${sample}.%j.out
#SBATCH -e $DEBUG/aligned_UTR_${sample}.%j.err
#SBATCH -t 10-00:00:00
#SBATCH -c 4
#SBATCH --mem-per-cpu=40G

set -eo pipefail   

module load python-anaconda3/2019.10
source /home/cpp4747/software/anaconda3/etc/profile.d/conda.sh
conda activate /home/cpp4747/software/anaconda3/envs/py39

echo "Running sample: \$sample"

find "\$sample_dir" -type f -name "*_utr.fasta" | while read utr_file; do
    fname=\$(basename "\$utr_file" .fasta)
    dedup_fa="\$dedup_dir/\$sample/\${fname}_dedup.fasta"

    if [[ ! -s "\$dedup_fa" ]]; then
        cp "\$utr_file" "\$dedup_fa"
    fi
done


for dedup_fa in "\$dedup_dir/\$sample"/*_dedup.fasta; do
    base=\$(basename "\$dedup_fa" _dedup.fasta)
	gene=\${base%%_*}
    aligned_fa="\$aligned_output/\$sample/\${base}_highscore_aligned_seqs.fa"
    unaligned_fa="\$aligned_output/\$sample/\${base}_unaligned.fa"

    echo "Aligning \$base"

    python3 "\$PY_SCRIPT" \
    "\$dedup_fa" \
    "\$GENE_REF_FA" \
    "\$aligned_fa" \
    "\$unaligned_fa" \
    --gene "\$gene" \
    --mismatch-rate 0.0 \
    --processes 4 \
    --batch-size 1000

done

echo "Finished sample: \$sample"
EOF

done
