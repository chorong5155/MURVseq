#!/bin/bash
topdir=$(pwd)
account="b1042"
queue="genomics"
DEBUG="$topdir/debug"
MAP="$topdir/gene_15mer_map.tsv"
refDir="/projects/b1154/referenceGenomes"
genomeID="hg38+VacV"
FASTA="${refDir}/${genomeID}/chrVacV_CDS15_with_phase.fa"
BASE="$topdir/polyA_output"

if [[ ! -f "$MAP" ]]; then
  awk '
  /^>/ {
    gsub(/^>/, "", $1)
    split($1, parts, "::")            
    split(parts[1], ginfo, "_")        
    gene = ginfo[1] "_" ginfo[2]       
    getline
    print gene "\t" $0
  }' "$FASTA" > "$MAP"
fi


declare -A SEQ15
while IFS=$'\t' read -r gene_phase seq; do
  SEQ15["$gene_phase"]="$seq"
done < "$MAP"

SEQ15_INIT=""
for key in "${!SEQ15[@]}"; do
  escaped_key=$(printf '%q' "$key")
  escaped_val=$(printf '%q' "${SEQ15[$key]}")
  SEQ15_INIT+="SEQ15[$escaped_key]=$escaped_val"$'\n'
done

for full_path in "$BASE"/*; do
  [[ -d "$full_path" ]] || continue
  sample=$(basename "$full_path")
  echo "Submitting job for sample: $sample"

  FASTQ_DIR="$full_path/polyA_matched_reads"
  OUTDIR="$full_path/polyA_length_per_gene"
  mkdir -p "$OUTDIR"

  sbatch <<EOF
#!/bin/bash -l
#SBATCH -A $account
#SBATCH -p $queue
#SBATCH -o $DEBUG/${sample}.out
#SBATCH -e $DEBUG/${sample}.err
#SBATCH -t 48:00:00
#SBATCH -c 2
#SBATCH --mem-per-cpu=40G
#SBATCH -J polyA_${sample}

FASTQ_DIR="$FASTQ_DIR"
OUTDIR="$OUTDIR"

declare -A SEQ15
$SEQ15_INIT

for gene_phase in "\${!SEQ15[@]}"; do
  seq15="\${SEQ15[\$gene_phase]}"
  gene="\${gene_phase%_*}"
  phase="\${gene_phase#*_}"

  fq="\${FASTQ_DIR}/\${gene}_\${phase}.fastq"
  out="\${OUTDIR}/\${gene}_\${phase}_polyA_lengths.txt"
  tmp="\${OUTDIR}/\${gene}_\${phase}.seqs"

  if [[ ! -f "\$fq" ]]; then
    echo "  Missing FASTQ for \$gene_phase: \$fq" >&2
    continue
  fi

  echo "  \$gene (\$phase): extracting polyA+\${#seq15}-mer…"

  awk 'NR%4==2' "\$fq" | grep -E -o "ATAT[A]{5,}\$seq15" > "\$tmp"

  total_length=0
  read_count=0
  declare -A length_counts
  while read -r match; do
    leader=\$(( \${#match} - 4 - \${#seq15} ))
    total_length=\$(( total_length + leader ))
    read_count=\$(( read_count + 1 ))
    length_counts[\$leader]=\$(( length_counts[\$leader] + 1 ))
  done < "\$tmp"

  if (( read_count > 0 )); then
    {
      echo "Poly-A leader length summary for \$gene (\$phase):"
      echo "Total reads: \$read_count"
      echo "Total poly-A leader length: \$total_length"
      echo
      echo "leader length counts:"
      for length in \$(printf "%s\n" "\${!length_counts[@]}" | sort -n); do
        printf "%4d %2d\n" "\${length_counts[\$length]}" "\$length"
      done
    } > "\$out"
    echo "    Wrote \$out"
  else
    echo "    No matches for \$gene (\$phase), removing empty output" >&2
    rm -f "\$out" "\$tmp"
  fi

  unset length_counts
done
EOF

done
