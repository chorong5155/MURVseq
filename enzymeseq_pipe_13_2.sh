#!/bin/bash
topdir=$(pwd)
account="p31754"
queue="long"

OUT_DIR="$topdir/classicslippage_TSS"
polyA="$topdir/polyA_output"
GENOME="/projects/b1154/referenceGenomes/hg38+VacV/hg38+VacV"
BED_FILE="/projects/b1154/referenceGenomes/hg38+VacV/VacV_150utr_flanks.bed"
BED_FILE_CDS="/projects/b1154/referenceGenomes/hg38+VacV/chrVacV_CDS15_with_phase.bed"
GENOME_FA="/projects/b1154/referenceGenomes/hg38+VacV/hg38+VacV.fa"

module load bedtools bowtie2 samtools gawk
mkdir -p "$OUT_DIR" "$topdir/debug"

for sample_dir in "$polyA"/*/; do
    sample=$(basename "$sample_dir")
    sample_out="$OUT_DIR/$sample"
    mkdir -p "$sample_out"

sbatch --export=ALL,SAMPLE="$sample",SAMPLE_DIR="$sample_dir",GENOME="$GENOME",BED_FILE="$BED_FILE",BED_FILE_CDS="$BED_FILE_CDS",GENOME_FA="$GENOME_FA",SAMPLE_OUT="$sample_out" <<EOF
#!/bin/bash -l
#SBATCH -A $account
#SBATCH -p $queue
#SBATCH -o $topdir/debug/TSS_${sample}.out
#SBATCH -e $topdir/debug/TSS_${sample}.err
#SBATCH -t 24:00:00
#SBATCH -c 2
#SBATCH --mem-per-cpu=40G
#SBATCH -J TSS_${sample}

module load bedtools samtools gawk bowtie2

EMPTY_REPORT="\$SAMPLE_OUT/empty_trimmed_counts.tsv"

for FASTA in "\$SAMPLE_DIR"/polyA_matched_reads/*_utr.fasta; do
    base=\$(basename "\$FASTA" .fasta)
    PREFIX="\$SAMPLE_OUT/\$base"

    # Trim polyA & assign empty reads to ORF start
awk -v BED="\$BED_FILE_CDS" -v PREFIX="\$PREFIX" -v gene="\$base" '
    BEGIN { OFS="\t"; seq=""; header="" }

    # When we hit a FASTA header line
    /^>/ {
        # Process previous read (if exists)
        if (header != "" && seq ~ /^A+$/) {
            readname = substr(header, 2)
            split(gene, arr, "_")
            gene_name = arr[1]

            # Load BED entry once
            if (!(gene_name in bed)) {
                while ((getline line < BED) > 0) {
                    split(line, b, "\t")
                    if (index(b[4], gene_name) > 0) {
                        bed[gene_name] = b[1] "\t" b[2] "\t" b[3] "\t" b[6]
                        break
                    }
                }
                close(BED)
            }

            # If gene exists in BED → write fallback TSS
            if (gene_name in bed) {
                split(bed[gene_name], bb, "\t")
                chr = bb[1]; start = bb[2]; end = bb[3]; strand = bb[4]

                if (strand == "+") tss = start
                else               tss = end

                bed_start = tss - 1
                if (bed_start < 0) bed_start = 0

                print chr, bed_start, tss, readname, 0, strand >> PREFIX "_ORF_fallback.bed"
            }
        }
        # Prepare for next record
        header = \$0
        seq = ""
        next
    }

    # Sequence lines
    /^[^>]/ { seq = seq \$0 }
END {
    # last read
    if(header!="") {
        if(seq !~ /^A+$/) print header "\n" seq >> PREFIX "_trimmed.fasta"
        else {
            # same polyA fallback as above
            readname = substr(header,2)
            split(gene, arr, "_")
            gene_name = arr[1]
            if (!(gene_name in bed)) {
                while ((getline line < BED) > 0) {
                    split(line,b,"\t")
                    if (index(b[4],gene_name)>0) { bed[gene_name]=b[1]"\t"b[2]"\t"b[3]"\t"b[6]; break }
                }
                close(BED)
            }
            if (gene_name in bed) {
                split(bed[gene_name], bb, "\t")
                chr=bb[1]; start=bb[2]; end=bb[3]; strand=bb[4]
                tss = (strand=="+")? start : end
                bed_start = tss-1; if(bed_start<0) bed_start=0
                print chr, bed_start, tss, readname, 0, strand >> PREFIX "_ORF_fallback.bed"
            }
        }
    }
}
' "\$FASTA"


    # Map non-empty reads if exist
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
                c=substr(cig,RSTART,RLENGTH);
                gsub(/[A-Z=]/,"",c);
                op=substr(cig,RSTART+length(c),1);
                n=c+0;
                if(op ~ /[MDN=X]/) refLen += n;
                cig=substr(cig,RSTART+RLENGTH);
            }
            if(flag % 32 >= 16) { strand="-"; tss=pos+refLen-1; }
            else { strand="+"; tss=pos; }
            start=tss-1; if(start<0) start=0;
            print chr "\t" start "\t" tss "\t" read "\t0\t" strand
        }' > "\${PREFIX}_TSS.bed"
    else
        echo "No non-empty reads for \$base, skipping mapping"
        touch "\${PREFIX}_TSS.bed"
    fi

    # Combine mapped TSS and ORF fallback
    COMBINED_BED="\${PREFIX}_TSS_combined.bed"
     #cut -f1-3 cat "\${PREFIX}_TSS.bed" "\${PREFIX}_ORF_fallback.bed" | sort -u > "\${PREFIX}_TSS_unique.bed"
 cat "\${PREFIX}_TSS.bed" "\${PREFIX}_ORF_fallback.bed" \
    | awk -v OFS="\t" '
        {
            chr=\$1; start=\$2; end=\$3; strand=\$6;
            key = chr "\t" start "\t" end "\t" strand;
            counts[key]++
        }
        END {
            for (k in counts) {
                # explode key back into fields
                split(k, a, "\t");
                chr=a[1]; start=a[2]; end=a[3]; strand=a[4];
                print chr, start, end, counts[k], strand
            }
        }
    ' \
    | sort -k1,1 -k2,2n \
    > "\${PREFIX}_TSS_counts.bed"
 


awk -v OFS="\t" '
{
    # TSS = end coordinate (3rd column)
    tss = \$3
    chr = \$1
    strand = \$5
    count = \$4

    key = chr "\t" tss "\t" strand

    counts[key] = count
}
END {
    for (k in counts) {
        print counts[k], k
    }
}
' "\${PREFIX}_TSS_counts.bed" \
| sort -nrk1,1 | head -3 \
| awk -v OFS="\t" '{print \$2, \$3, \$4, \$1}' \
> "\${PREFIX}_dominant_TSS.txt"



done
EOF
done
