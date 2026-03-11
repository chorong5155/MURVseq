#!/bin/bash
topdir=$(pwd)
account="p31754"
queue="long"

OUT_DIR="$topdir/natural_TSS"
natural="$topdir/un_natural_output"
GENOME="/projects/b1154/referenceGenomes/hg38+VacV/hg38+VacV"
BED_FILE="/projects/b1154/referenceGenomes/hg38+VacV/VacV_150utr_flanks.bed"
GENOME_FA="/projects/b1154/referenceGenomes/hg38+VacV/hg38+VacV.fa"

module load bedtools bowtie2 samtools
mkdir -p "$OUT_DIR" "$topdir/debug"

for sample_dir in "$natural"/*/; do
    sample=$(basename "$sample_dir")
    sample_out="$OUT_DIR/$sample"
    mkdir -p "$sample_out"

    sbatch --export=ALL,SAMPLE="$sample",SAMPLE_DIR="$sample_dir",OUT_DIR="$sample_out",GENOME="$GENOME",BED_FILE="$BED_FILE",GENOME_FA="$GENOME_FA" <<EOF
#!/bin/bash -l
#SBATCH -A $account
#SBATCH -p $queue
#SBATCH -o $topdir/debug/TSS_${sample}.out
#SBATCH -e $topdir/debug/TSS_${sample}.err
#SBATCH -t 24:00:00
#SBATCH -c 2
#SBATCH --mem-per-cpu=40G
#SBATCH -J TSS_${sample}

module load bedtools samtools minimap2

EMPTY_REPORT="\$OUT_DIR/empty_trimmed_counts.tsv"

for FASTA in "\$SAMPLE_DIR"/naturalUTRs/*_naturalUTRs.fasta; do
    base=\$(basename "\$FASTA" .fasta)
    PREFIX="\$OUT_DIR/\$base"

   
    bowtie2 --local -L 14 -f -x "\$GENOME" -U "\${FASTA}" -S "\${PREFIX}.sam"

    samtools view -bS "\${PREFIX}.sam" | samtools sort -o "\${PREFIX}.sorted.bam"
    samtools index "\${PREFIX}.sorted.bam"
    samtools view "\${PREFIX}.sorted.bam" |
    awk '
        # Skip unmapped reads
        (\$2 % 8) >= 4 { next }

        {
            chr=\$3; pos=\$4; read=\$1; flag=\$2; cigar=\$6

            # compute reference-consuming length
            refLen=0
            cig=cigar
            while (match(cig, /[0-9]+[MIDNSHP=X]/)) {
                c=substr(cig, RSTART, RLENGTH)
                gsub(/[A-Z=]/, "", c)
                op=substr(cig, RSTART+length(c), 1)
                n=c+0

                if (op=="M" || op=="D" || op=="N" || op=="=" || op=="X")
                    refLen += n

                cig=substr(cig, RSTART+RLENGTH)
            }

            # strand logic (poxvirus)
            if (flag % 32 >= 16) {
                strand = "-"
                tss = pos + refLen - 1
            } else {
                strand = "+"
                tss = pos
            }

            start = tss - 1
            if (start < 0) start = 0

            print chr "\t" start "\t" tss "\t" read "\t.\t" strand
        }
    ' > "\${PREFIX}_TSS.bed" 

awk -v OFS="\t" '{ key=\$1"\t"\$2"\t"\$3; counts[key]++ }
END {
    for (k in counts) print k, counts[k]
}' "\${PREFIX}_TSS.bed" \
| sort -k1,1 -k2,2n > "\${PREFIX}_TSS_unique.bed"


    bedtools intersect -a "\${PREFIX}_TSS.bed" -b "\$BED_FILE" -wa -wb \
        > "\${PREFIX}_TSS_annotated.bed"

    #cut -f1-3 "\${PREFIX}_TSS.bed" | sort -u > "\${PREFIX}_TSS_unique.bed"

awk -v OFS="\t" '
{
    key=\$1"\t"\$3"\t"\$6
    counts[key]++
}
END {
    for (k in counts) {
        print counts[k], k
    }
}
' "\${PREFIX}_TSS.bed" \
| sort -nrk1,1 | head -3 \
| awk -v OFS="\t" '{print \$2, \$3, \$4, \$1}' \
> "\${PREFIX}_dominant_TSS.txt"




    bedtools getfasta \
        -fi "\${GENOME_FA}" \
        -bed "\${PREFIX}_TSS.bed" \
        -s \
        -name \
        -fo "\${PREFIX}_TSS_seq.fa"

    awk -v OFS="\t" -v stage="early" '{
        TSS_ID="TSS_"NR
        print \$1, \$2, \$3, TSS_ID, ".", \$6, \$7, \$8, \$9, \$10, stage
    }' "\${PREFIX}_TSS_annotated.bed" > "\${PREFIX}_TSS_final.tsv"


done
EOF
done
