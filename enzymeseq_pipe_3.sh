#!/bin/bash
topdir=$(pwd)
account="b1042"
queue="genomics"
outdir="$topdir/debug"
outputdir="$topdir/polyA_output_2"
refDir="/projects/b1154/referenceGenomes"
genomeID="hg38+VacV"
mkdir -p "$outdir" "$outputdir"
output_root="$topdir/output"
emailaddress="${USER}@e.northwestern.edu"
CDS15="${refDir}/${genomeID}/chrVacV_CDS15_with_phase.fa"
if [[ ! -s "$CDS15" ]]; then
    echo "ERROR: File not found or empty: $CDS15" >&2
    exit 1
fi

echo "Searching FASTQ files under: $output_root"

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

    MATCHDIR="$outputdir/${sample}/polyA_matched_reads"
    mkdir -p "$MATCHDIR"


    sbatch <<EOF
#!/bin/bash -l
#SBATCH -A $account
#SBATCH -p $queue
#SBATCH -o $outdir/polyA_${sample}.out
#SBATCH -e $outdir/polyA_${sample}.err
#SBATCH -t 48:00:00
#SBATCH -c 2
#SBATCH --mem-per-cpu=40G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=$emailaddress
#SBATCH -J polyA_${sample}

echo "Starting job for $sample at: \$(date)"
mkdir -p "$MATCHDIR"

# Step 1: extract gene info into a temporary gene list
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
        regex = "ATATAA*A*A{3,}" seq15;
        filename = gene "_" phase;
        print gene, phase, regex, filename, seq15;
        nextline = 0;
    }' "$CDS15" > "$MATCHDIR/gene_list.txt"

# Step 2: loop through gene list and process FASTQ reads
while read -r gene phase regex filename seq15; do
    outfile="$MATCHDIR/\${filename}.fastq"
    utrout="$MATCHDIR/\${filename}_utr.fasta"

    # Skip if output already exists
    if [[ -s "\$outfile" ]]; then
        echo "Skipping \$outfile — already exists"
        continue
    fi

    for fastq in "$r1" "$r2"; do
        if [[ ! -f "\$fastq" ]]; then
            echo "WARNING: File not found: \$fastq"
            continue
        fi

        awk -v pattern="\$regex" 'BEGIN {OFS="\n"}
        NR % 4 == 1 {h=\$0}
        NR % 4 == 2 {s=\$0}
        NR % 4 == 3 {p=\$0}
        NR % 4 == 0 {
            q=\$0
            if (s ~ pattern) {
                print h,s,p,q
            }
        }' "\$fastq" >> "\$outfile"
    done 
    
    # Extract UTR sequences from matched reads
    awk -v seq15="\$seq15" 'BEGIN {OFS="\n"}
    NR % 4 == 1 {header=substr(\$0, 2)}
    NR % 4 == 2 {
        read=\$0
        if (match(read, /(ATATAA+)(A*)/, arr)) {
            atat_start=RSTART
            atat_end=RSTART + 3
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

done < "$MATCHDIR/gene_list.txt"

echo "Finished job for $sample at: \$(date)"
EOF

done
