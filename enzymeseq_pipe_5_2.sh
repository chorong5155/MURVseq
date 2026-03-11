#!/bin/bash
topdir=$(pwd)
account="b1042"
queue="long"
refDir="/projects/b1154/referenceGenomes"
genomeID="hg38+VacV"
CDS15="${refDir}/${genomeID}/chrVacV_CDS15_with_phase.fa"
output_root="$topdir/output"
fastq_dir="${topdir}/fastq/unzipped_files"
total_output="${topdir}/total_output"
polya_output="${topdir}/polyA_output"
DEBUG="${topdir}/debug"

mkdir -p "$total_output" "$polya_output" "$DEBUG"
CDS15="${refDir}/${genomeID}/chrVacV_CDS15_with_phase.fa"
if [[ ! -s "$CDS15" ]]; then
    echo "ERROR: File not found or empty: $CDS15" >&2
    exit 1
fi

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

    echo "Sample: $sample"
    echo "R1: $r1"
    echo "R2: $r2"

    sample_out="$total_output/$sample"
    mkdir -p "$sample_out"

sbatch <<EOF
#!/bin/bash -l
#SBATCH -A b1042
#SBATCH --partition=genomicslong
#SBATCH -o $DEBUG/total_${sample}.out
#SBATCH -e $DEBUG/total_${sample}.err
#SBATCH -t 10-00:00:00
#SBATCH -c 2
#SBATCH --mem-per-cpu=40G
#SBATCH -J total_${sample}

echo "Starting job for $sample at: \$(date)"
mkdir -p "$sample_out"

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
    regex = "ATAT[ACGTN]{1,130}" seq15;
    filename = gene "_" phase;
    print gene, phase, regex, filename, seq15;
    nextline = 0;
}' "$CDS15" > "$sample_out/gene_list.txt"

 while read -r gene phase regex filename seq15; do
    
    outfile="${sample_out}/\${filename}.fastq"
	utrout="${sample_out}/\${filename}_utr.fasta"
	
	if [[ ! -s "\$outfile" ]]; then
       echo "Creating \$outfile"
    for fq in "$r1" "$r2"; do
        [[ -f "\$fq" ]] || continue
        awk -v pattern="\$regex" 'BEGIN {OFS="\n"}
        BEGIN {OFS="\\n"}
        NR % 4 == 1 {h=\$0}
        NR % 4 == 2 {s=\$0}
        NR % 4 == 3 {p=\$0}
        NR % 4 == 0 {
            q=\$0
            if (s ~ pattern) {
                print h,s,p,q
            }
        }' "\$fq" >> "\$outfile"
        done
        else
    echo "Skipping \$outfile — already exists"
fi

  awk -v seq15="\$seq15" 'BEGIN {OFS="\n"}
    NR % 4 == 1 {header=substr(\$0, 2)}
    NR % 4 == 2 {
        read=\$0
        if (match(read, /ATAT[ACGTN]{1,130}/)) {
            atat_start = RSTART;
             atat_end=RSTART + 3;
             } else {
            atat_start=0
            atat_end=0
        }
        cds_start=index(read, seq15)
        print header, "ATATA end:", atat_end, "CDS start:", cds_start > "/dev/stderr"
        if (atat_start > 0 && cds_start > 0 && cds_start > atat_end) {
            utr=substr(read, atat_end+1, cds_start-atat_end-1)
            print ">" header, utr
        }
    }' "\$outfile" > "\$utrout"
    
    done < "$sample_out/gene_list.txt"

echo "Finished job for $sample at: \$(date)"
EOF

done
