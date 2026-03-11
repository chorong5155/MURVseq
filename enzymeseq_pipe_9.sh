#!/bin/bash
set -euo pipefail
topdir=$(pwd)
account="b1042"
queue="genomics"
fastq_dir="${topdir}/output"
DEBUG="${topdir}/debug"

for sampleDir in "$fastq_dir"/*/; do
    sample=$(basename "$sampleDir")

    r1="$sampleDir/${sample}_1P_czid_dedup.fastq"
    r2="$sampleDir/${sample}_2P_czid_dedup.fastq"

    if [[ ! -f "$r1" ]]; then
        echo "Missing R1 for sample: $sample"
        continue
    fi
    if [[ ! -f "$r2" ]]; then
        echo "Missing R2 for sample: $sample"
        continue
    fi

 mkdir -p "$DEBUG"
  
  NATURAL_DIR="$topdir/un_natural_output/${sample}/naturalUTRs"
  UNNATURAL_DIR="$topdir/un_natural_output/${sample}/unnaturalUTRs"
  POLYA_DIR="$topdir/polyA_output/${sample}/polyA_matched_reads"
  MICRO_DIR="$topdir/micropolyA_dedup/${sample}"
  TOTAL_DIR="$topdir/total_output/${sample}"
  OUTPUT_DIR="$topdir/proportion_output/${sample}"
  mkdir -p "$OUTPUT_DIR"
  
  sbatch <<EOF
#!/bin/bash -l
#SBATCH -A $account
#SBATCH -p $queue
#SBATCH -o $topdir/debug/prop_${sample}.out
#SBATCH -e $topdir/debug/prop_${sample}.err
#SBATCH -t 48:00:00
#SBATCH -c 2
#SBATCH --mem-per-cpu=10G
#SBATCH -J prop${sample}

# Debug prints to check paths
echo "Sample: $sample"
echo "NATURAL_DIR=$NATURAL_DIR"
echo "UNNATURAL_DIR=$UNNATURAL_DIR"
echo "POLYA_DIR=$POLYA_DIR"
echo "TOTAL_DIR=$TOTAL_DIR"
echo "MICRO_DIR=$MICRO_DIR"
echo "OUTPUT_DIR=$OUTPUT_DIR"

count_headers_per_gene() {
  local dir="\$1"
  local outfile="\$2"

  if [ ! -d "\$dir" ]; then
    echo "Warning: directory \$dir does not exist. Skipping."
    echo "" > "\$outfile"
    return
  fi

  declare -A gene_counts

  for file in "\$dir"/*.fasta; do
    [ -e "\$file" ] || continue
    echo "Processing \$file"
    filename=\$(basename "\$file")
    gene=\$(echo "\$filename" | cut -d'_' -f1)
    phase=\$(echo "\$filename" | cut -d'_' -f2)
    key="\${gene}_\${phase}"
    count=\$(grep -c "^>" "\$file")
    ((gene_counts["\$key"] += count))
  done

  {
    for key in "\${!gene_counts[@]}"; do
      echo -e "\$key\t\${gene_counts[\$key]}"
    done | sort
  } > "\$outfile"

  echo "Counted headers in \$dir -> \$outfile"
}

# Count fasta headers using fully expanded paths
count_headers_per_gene "$NATURAL_DIR" "$OUTPUT_DIR/counts_natural.txt"
count_headers_per_gene "$UNNATURAL_DIR" "$OUTPUT_DIR/counts_unnatural.txt"
count_headers_per_gene "$POLYA_DIR" "$OUTPUT_DIR/counts_polya.txt"
count_headers_per_gene "$TOTAL_DIR" "$OUTPUT_DIR/counts_total.txt"
count_headers_per_gene "$MICRO_DIR" "$OUTPUT_DIR/counts_micropolyA.txt"

echo "Counting done!"

for f in counts_natural counts_unnatural counts_polya counts_micropolyA counts_total; do  # include micropolyA
  awk '{print \$1, \$2}' "$OUTPUT_DIR/\${f}.txt" > "$OUTPUT_DIR/\${f}.tmp"
done


cat "$OUTPUT_DIR"/*.tmp | cut -d' ' -f1 | sort | uniq > "$OUTPUT_DIR/all_genes.txt"

# Header
echo -e "Gene_Timepoint\tNatural\tUnnatural\tpolyA\tmicropolyA\ttotal\tprop_polyA\tprop_micropolyA\tprop_natural" > "$OUTPUT_DIR/combined_wide.txt"


while read gene; do
  n=\$(grep -w "^\$gene" "$OUTPUT_DIR/counts_natural.tmp" | cut -d' ' -f2)
  u=\$(grep -w "^\$gene" "$OUTPUT_DIR/counts_unnatural.tmp" | cut -d' ' -f2)
  p=\$(grep -w "^\$gene" "$OUTPUT_DIR/counts_polya.tmp" | cut -d' ' -f2)
  m=\$(grep -w "^\$gene" "$OUTPUT_DIR/counts_micropolyA.tmp"| cut -d' ' -f2) 
  t=\$(grep -w "^\$gene" "$OUTPUT_DIR/counts_total.tmp" | cut -d' ' -f2)

  n=\${n:-0}
  u=\${u:-0}
  p=\${p:-0}
  m=\${m:-0}
  t=\${t:-0}

 u=\$(( u - (p + m) ))
  if [ "\$u" -lt 0 ]; then
    u=0   # clamp to 0 if subtraction goes negative
  fi
  
if [ "\$t" -gt 0 ]; then
    prop_n=\$(awk -v n="\$n" -v t="\$t" 'BEGIN{printf "%.4f", n/t}')
    prop_p=\$(awk -v p="\$p" -v t="\$t" 'BEGIN{printf "%.4f", p/t}')
    prop_m=\$(awk -v m="\$m" -v t="\$t" 'BEGIN{printf "%.4f", m/t}')
else
    prop_n=0
    prop_p=0
    prop_m=0
fi
  echo -e "\$gene\t\$n\t\$u\t\$p\t\$m\t\$t\t\$prop_p\t\$prop_m\t\$prop_n"
done < "$OUTPUT_DIR/all_genes.txt" >> "$OUTPUT_DIR/combined_wide.txt"

# Clean up
rm "$OUTPUT_DIR"/*.tmp
rm "$OUTPUT_DIR/all_genes.txt"

EOF

done
