#!/bin/bash
topdir=$(pwd)
account="p31754"
queue="long"

OUT_DIR="$topdir/slippage_TSS"
micropolyA="$topdir/micropolyA_dedup"

GENOME="/projects/b1154/referenceGenomes/hg38+VacV/hg38+VacV"
BED_FILE="/projects/b1154/referenceGenomes/hg38+VacV/VacV_150utr_flanks.bed"
BED_FILE_CDS="/projects/b1154/referenceGenomes/hg38+VacV/chrVacV_CDS15_with_phase.bed"
GENOME_FA="/projects/b1154/referenceGenomes/hg38+VacV/hg38+VacV.fa"

module load bedtools bowtie2 samtools
mkdir -p "$OUT_DIR" "$topdir/debug"

for sample_dir in "$micropolyA"/*/; do
    sample=$(basename "$sample_dir")
    sample_out="$OUT_DIR/$sample"
    mkdir -p "$sample_out"

sbatch --export=ALL,SAMPLE="$sample",SAMPLE_DIR="$sample_dir",SAMPLE_OUT="$sample_out",OUT_DIR="$sample_out",GENOME="$GENOME",BED_FILE="$BED_FILE",BED_FILE_CDS="$BED_FILE_CDS",GENOME_FA="$GENOME_FA" <<EOF
#!/bin/bash -l
#SBATCH -A $account
#SBATCH -p $queue
#SBATCH -o $topdir/debug/TSSmicro_${sample}.out
#SBATCH -e $topdir/debug/TSSmicro_${sample}.err
#SBATCH -t 24:00:00
#SBATCH -c 2
#SBATCH --mem-per-cpu=40G
#SBATCH -J TSSmicro_${sample}

module load bedtools samtools minimap2 bowtie2

for FASTA in "\$SAMPLE_DIR"/*_dedup.fasta; do
    base=\$(basename "\$FASTA" .fasta)
    PREFIX="\$OUT_DIR/\$base"

    awk -v BED="\$BED_FILE_CDS" -v PREFIX="\$PREFIX" -v gene="\$base" '
    BEGIN {
        OFS="\t"; seq=""; header="";
        while ((getline line < BED) > 0) {
            split(line, b, "\t")
            bed[b[4]] = b[1] "\t" b[2] "\t" b[3] "\t" b[6]
        }
        close(BED)
    }

    /^>/ {
        if (header != "") {
            readname = substr(header, 2)
            split(gene, arr, "_")
            gene_name = arr[1]

            if (seq ~ /^A+\$/ && (gene_name in bed)) {
                split(bed[gene_name], bb, "\t")
                chr = bb[1]; start = bb[2]; end = bb[3]; strand = bb[4]
                tss = (strand == "+") ? start : end
                bed_start = (tss > 0 ? tss - 1 : 0)
                print chr, bed_start, tss, readname, 0, strand >> PREFIX "_ORF_fallback.bed"
            } else {
                print header "\n" seq >> PREFIX "_trimmed.fasta"
            }
        }
        header = \$0; seq=""; next
    }

    /^[^>]/ { seq = seq \$0 }

    END {
        if (header != "") {
            readname = substr(header, 2)
            split(gene, arr, "_")
            gene_name = arr[1]
            if (seq ~ /^A+\$/ && (gene_name in bed)) {
                split(bed[gene_name], bb, "\t")
                chr = bb[1]; start = bb[2]; end = bb[3]; strand = bb[4]
                tss = (strand == "+") ? start : end
                bed_start = (tss > 0 ? tss - 1 : 0)
                print chr, bed_start, tss, readname, 0, strand >> PREFIX "_ORF_fallback.bed"
            } else {
                print header "\n" seq >> PREFIX "_trimmed.fasta"
            }
        }
    }
    ' "\$FASTA"

    if [ -s "\${PREFIX}_trimmed.fasta" ]; then
        bowtie2 --local -L 14 -f -x "\$GENOME" -U "\${PREFIX}_trimmed.fasta" -S "\${PREFIX}.sam"
        samtools view -bS "\${PREFIX}.sam" | samtools sort -o "\${PREFIX}.sorted.bam"
        samtools index "\${PREFIX}.sorted.bam"

        samtools view "\${PREFIX}.sorted.bam" | awk '
        (\$2 % 8) >= 4 { next }
        {
            chr=\$3; pos=\$4; read=\$1; flag=\$2; cigar=\$6;
            refLen=0; cig=cigar;
            while(match(cig, /[0-9]+[MIDNSHP=X]/)) {
                c=substr(cig,RSTART,RLENGTH)
                gsub(/[A-Z=]/,"",c)
                op=substr(cig,RSTART+length(c),1)
                n=c+0
                if(op ~ /[MDN=X]/) refLen += n
                cig=substr(cig,RSTART+RLENGTH)
            }
            if(flag % 32 >= 16) { strand="-"; tss=pos+refLen-1 }
            else { strand="+"; tss=pos }
            start = (tss>0 ? tss-1 : 0)
            print chr "\t" start "\t" tss "\t" read "\t0\t" strand
        }' > "\${PREFIX}_TSS.bed"

    else
        echo "No non-empty reads for \$base"
        touch "\${PREFIX}_TSS.bed"
    fi

    awk -v OFS="\t" '
    {
        key = \$1"\t"\$3"\t"\$6
        counts[key]++
    }
    END {
        for (k in counts) {
            split(k, a, "\t")
            start = a[2]-1; if(start<0) start=0
            print a[1], start, a[2], counts[k], 0
        }
    }' "\${PREFIX}_TSS.bed" > "\${PREFIX}_TSS_counts.bed"

    awk -v OFS="\t" '
    {
        key=\$1"\t"\$3"\t"\$6
        counts[key]++
    }
    END {
        for (k in counts)
            print counts[k], k
    }' "\${PREFIX}_TSS.bed" \
    | sort -nrk1,1 | head -3 \
     | awk -v OFS="\t" '{split(\$0,a," "); print a[2], a[3], a[4], a[1]}' \
    > "\${PREFIX}_dominant_TSS.txt"

    cat "\${PREFIX}_TSS.bed" "\${PREFIX}_ORF_fallback.bed" \
    | cut -f1-3 | sort -u \
    > "\${PREFIX}_TSS_unique.bed"

    cat "\${PREFIX}_TSS.bed" "\${PREFIX}_ORF_fallback.bed" \
    | cut -f1-3 | sort -u \
    > "\${PREFIX}_TSS_unique.bed"



done
EOF
done
