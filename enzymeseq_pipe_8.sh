#!/bin/bash
topdir=$(pwd)
account="b1042"
queue="genomics"
GENOME="/projects/b1154/referenceGenomes/hg38+VacV"
CDS15="$GENOME/chrVacV_CDS15_with_phase.fa"
DEBUG="${topdir}/debug"
output_root="$topdir/output"
total_dir="${topdir}/total_output"

for sampleDir in "$output_root"/*/; do
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

 MATCHDIR="${topdir}/micropolyA_output/${sample}"
 mkdir -p "$MATCHDIR" "$DEBUG"
 
sbatch <<EOF
#!/bin/bash -l
#SBATCH -A $account
#SBATCH -p $queue
#SBATCH -o $DEBUG/micro_${sample}.out
#SBATCH -e $DEBUG/micro_${sample}.err
#SBATCH -t 48:00:00
#SBATCH -c 2
#SBATCH --mem-per-cpu=40G
#SBATCH -J micropolyA${sample}


echo "Starting job for $sample at: \$(date)"

awk '/^>/ {
        sub(/^>/, "", \$0);
        split(\$0, parts, "::");
        split(parts[1], g, "_");
        gene = g[1];
        phase = g[2];
        nextline = 1;
        next;
    }
nextline {
    seq15 = \$0;
    regex = "ATATAA{3,}[ACGT]{1,}" seq15;
    filename = gene "_" phase;
    print gene, phase, regex, filename, seq15;
    nextline = 0;
}' "$CDS15" > "$MATCHDIR/gene_list.txt"
	
	while read gene phase regex filename seq15; do
    polyA_outfile="$MATCHDIR/\${filename}.fastq"
    utrout="$MATCHDIR/\${filename}_utr.fasta"
    echo "Searching for \$gene (\$phase) with pattern: \$regex"

   if [[ -s "\$polyA_outfile" ]]; then
    echo "  skipping, \$polyA_outfile already exists...."
	else
    for fastq in "$r1" "$r2"; do
        awk -v r="\$regex" 'BEGIN{OFS="\n"} 
            NR%4==1 {h=\$0} 
            NR%4==2 {s=\$0} 
            NR%4==3 {p=\$0} 
            NR%4==0 {q=\$0; if (s ~ r) print h,s,p,q}' "\$fastq" >> "\$polyA_outfile"
     done
    
fi

    if [[ ! -s "\$utrout" ]]; then
        echo "  Extracting UTRs into \$utrout..."
        awk -v seq15="\$seq15" 'BEGIN {OFS="\n"} 
            NR % 4 == 1 {header=substr(\$0, 2)} 
            NR % 4 == 2 { 
                read=\$0; 
                if (match(read, /ATATAA{3,}[ACGT]{1,}/)) {
                    atat_start = RSTART; 
                    atat_end = RSTART + 3; 
                } else { 
                    atat_start = 0; 
                    atat_end = 0; 
                } 
                cds_start = index(read, seq15); 
                print header, "ATATA end:", atat_end, "CDS start:", cds_start > "/dev/stderr"; 
                if (atat_start > 0 && cds_start > 0 && cds_start > atat_end) { 
                    utr = substr(read, atat_end + 1, cds_start - atat_end - 1); 
                    print ">" header, utr; 
                } 
            }' "\$polyA_outfile" > "\$utrout"
    else
        echo "  Skipping UTR extraction, \$utrout already exists."
    fi

done < "$MATCHDIR/gene_list.txt"

echo "Finished job for $sample at: \$(date)"
EOF
done