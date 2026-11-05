#!/bin/bash
topdir=$(pwd)
account="b1042"
queue="genomics"
BASE_DIR="/projects/b1042/WalshLab/cpp4747/19215FL-11/aligned_UTR"
CDS_FASTA="/projects/b1154/cpp4747/hg38+VacV/chrVacV_fullCDS.fa"
out_dir="${topdir}/uAUGtest"
DEBUG="${topdir}/debug"

mkdir -p "$out_dir" "$DEBUG"

echo "OUTER SCRIPT STARTED"
echo "topdir     = $topdir"
echo "BASE_DIR   = $BASE_DIR"
echo "CDS_FASTA  = $CDS_FASTA"
echo "out_dir    = $out_dir"
echo "DEBUG      = $DEBUG"

if [[ ! -d "$BASE_DIR" ]]; then
    echo "ERROR: BASE_DIR does not exist:"
    echo "$BASE_DIR"
    exit 1
fi

if [[ ! -f "$CDS_FASTA" ]]; then
    echo "ERROR: CDS FASTA does not exist:"
    echo "$CDS_FASTA"
    exit 1
fi

echo "BASE_DIR and CDS FASTA exist."

sample_count=0

for sample_dir in "$BASE_DIR"/*; do
    [[ -d "$sample_dir" ]] || continue

    sample=$(basename "$sample_dir")
    sample_count=$((sample_count + 1))
    sample_out="$out_dir/$sample"

    mkdir -p "$sample_out"
    fa_count=$(find "$sample_dir" -maxdepth 1 -type f -name "*.fa" | wc -l)

    echo
    echo "Sample: $sample"
    echo "FASTA files found: $fa_count"

    if [[ "$fa_count" -eq 0 ]]; then
        echo "WARNING: No FASTA files found. Skipping $sample."
        continue
    fi

sbatch --export=ALL,sample_dir="$sample_dir",sample="$sample",sample_out="$sample_out" <<EOF
#!/bin/bash -l
#SBATCH -A $account
#SBATCH -p $queue
#SBATCH -o $DEBUG/uAUG_all_${sample}.out
#SBATCH -e $DEBUG/uAUG_all_${sample}.err
#SBATCH -t 48:00:00
#SBATCH -c 2
#SBATCH --mem-per-cpu=40G
#SBATCH -J uAUG_all${sample}

STOP_CODONS=("TAA" "TAG" "TGA")

echo "SLURM JOB STARTED"
echo "Sample: \$sample"
echo "Date: \$(date)"

if [[ ! -d "\$sample_dir" ]]; then
    echo "ERROR: sample directory does not exist:"
    echo "\$sample_dir"
    exit 1
fi

if [[ ! -f "$CDS_FASTA" ]]; then
    echo "ERROR: CDS FASTA does not exist:"
    echo "$CDS_FASTA"
    exit 1
fi

mkdir -p "\$sample_out"

# Build CDS 15-mer lookup file (extracts first 15 bp of CDS)
cds_lookup="\$sample_out/CDS_lookup.tsv"
echo "Building CDS 15-mer lookup..."

awk '
    /^>/ {
        if (gene != "")
            print gene "\t" substr(seq, 1, 15)

        gene=\$1
        sub(/^>/, "", gene)
        sub(/::.*/, "", gene)
        seq=""
        next
    }
    {
        gsub(/[[:space:]\r]/, "", \$0)
        seq=seq \$0
    }
    END {
        if (gene != "")
            print gene "\t" substr(seq, 1, 15)
    }
' "$CDS_FASTA" > "\$cds_lookup"

# find CDS 15-mer sequence

get_cds() {
    local gene="\$1"

    awk -v g="\$gene" '
        \$1 == g {
            print \$2
            exit
        }
    ' "\$cds_lookup"
}

# find canonical ATG0-frame STOP in CDS 15-mer
# Returns 1-based position of the first STOP codon, or NA.

get_canonical_stop_pos() {
    local cds_seq="\$1"

    if [[ -z "\$cds_seq" ]]; then
        echo "NA"
        return
    fi

    local cds_len=\${#cds_seq}
    local j codon stop

    for ((j=0; j+2<\$cds_len; j+=3)); do

        codon=\${cds_seq:\$j:3}

        for stop in "\${STOP_CODONS[@]}"; do
            if [[ "\$codon" == "\$stop" ]]; then
                echo \$((j + 1))
                return
            fi
        done
    done

    echo "NA"
}

# classify one FASTA record

classify_sequence() {

    local header="\$1"
    local seq="\$2"
    local gene="\$3"
    local cds_seq="\$4"
    local canonical_stop_pos="\$5"
    local output_file="\$6"

    seq=\$(printf '%s' "\$seq" | tr '[:lower:]' '[:upper:]')

    # Defining ATG0 

    seq_len=\${#seq}
    ATG0_POS=\$((seq_len + 1))

    # Find first ATG in UTR = ATG1
    # grep position is 0-based, convert to 1-based

    ATG1_POS=\$(printf '%s' "\$seq" |
        grep -b -o "ATG" |
        head -n1 |
        cut -d: -f1)

    # Category 4: no uAUG

    if [[ -z "\$ATG1_POS" ]]; then
        echo -e "\$header\tNA\t\$ATG0_POS\tNA\tNA\tNA\tNA\tCategory4_No_uAUG" \
            >> "\$output_file"
        return
    fi

    ATG1_POS=\$((ATG1_POS + 1))

    # Determine ATG1 frame relative to ATG0 (Important for c1 and c3)
    diff=\$((ATG0_POS - ATG1_POS))
    frame=\$((diff % 3))

    # Scan from ATG1 toward ATG0 in ATG1 reading frame
    utr_subseq=\${seq:\$((ATG1_POS - 1))}
    subseq_len=\${#utr_subseq}

    stop_found=0

    if [[ \$subseq_len -ge 3 ]]; then

        for ((i=0; i+2<\$subseq_len; i+=3)); do

            codon=\${utr_subseq:\$i:3}

            for stop in "\${STOP_CODONS[@]}"; do

                if [[ "\$codon" == "\$stop" ]]; then
                    stop_found=1
                    break 2
                fi

            done
        done
    fi

    # Category 2: First STOP in the ATG1 reading frame occurs in the UTR before canonical ATG0.
    atg1_cds_stop_pos="NA"

    if [[ "\$stop_found" -eq 1 ]]; then
        category="Category2_Upstream_Stopped"

    # No STOP before ATG0.
    elif [[ "\$frame" -eq 0 ]]; then
        
        if [[ -n "\$cds_seq" && "\$canonical_stop_pos" != "NA" ]]; then
            atg1_cds_stop_pos="\$canonical_stop_pos"
            category="Category1_Productive_Upstream"
        else
            category="Category5"
        fi

    else

        if [[ -z "\$cds_seq" ]]; then
            category="Category5"
        else
            cds_stop_found=0
            cds_frame_offset=\$(((3 - frame) % 3))
            cds_len=\${#cds_seq}

            for ((i=\$cds_frame_offset; i+2<\$cds_len; i+=3)); do

                codon=\${cds_seq:\$i:3}

                for stop in "\${STOP_CODONS[@]}"; do

                    if [[ "\$codon" == "\$stop" ]]; then
                        cds_stop_found=1
                        # 1-based position within canonical CDS 15-mer.
                        atg1_cds_stop_pos=\$((i + 1))
                        break 2
                    fi

                done
            done

            if [[ "\$cds_stop_found" -eq 1 ]]; then
                category="Category3_Overlapping_OutOfFrame"
            else
                category="Category5"
            fi

        fi
    fi

    echo -e "\$header\t\$ATG1_POS\t\$ATG0_POS\t\$diff\t\$frame\t\${canonical_stop_pos:-NA}\t\${atg1_cds_stop_pos:-NA}\t\$category" \
        >> "\$output_file"
}

# PROCESS ALIGNED FASTA

aligned_count=0

for fa in "\$sample_dir"/*_highscore_aligned_seqs.fa; do

    [[ -f "\$fa" ]] || continue

    aligned_count=\$((aligned_count + 1))

    filename=\$(basename "\$fa")
    base=\${filename%_highscore_aligned_seqs.fa}

    # Gene comes from the FASTA filename
    gene=\${base%%_*}

    # Retrieve the matching canonical CDS 15-mer ONCE for this gene/file.
    cds_seq=\$(get_cds "\$gene")
    cds_seq=\$(printf '%s' "\$cds_seq" | tr '[:lower:]' '[:upper:]')
    canonical_stop_pos=\$(get_canonical_stop_pos "\$cds_seq")

    output_file="\$sample_out/\${base}_highscore_aligned_seqs_categorized.tsv"

    echo
    echo "Processing ALIGNED FASTA"
    echo "File: \$filename"
    echo "Gene: \$gene"
    echo "CDS 15-mer length: \${#cds_seq}"
    echo "Canonical CDS 15-mer STOP pos: \$canonical_stop_pos"
    echo "Output: \$output_file"

    echo -e "header\tatg1_pos\tatg0_pos\tdiff\tframe\tcanonical_cds_stop_pos\tatg1_frame_cds_stop_pos\tcategory" \
        > "\$output_file"

    awk '
        /^>/ {
            if (header != "")
                print header "\t" seq

            header=\$0
            seq=""
            next
        }
        {
            gsub(/[[:space:]\r]/, "", \$0)
            seq=seq \$0
        }
        END {
            if (header != "")
                print header "\t" seq
        }
    ' "\$fa" |
    while IFS=\$'\t' read -r header seq; do

        classify_sequence "\$header" "\$seq" "\$gene" "\$cds_seq" "\$canonical_stop_pos" "\$output_file"

    done

done

echo
echo "Aligned FASTA files processed: \$aligned_count"

# PROCESS UNALIGNED FASTA

unaligned_count=0

for fa in "\$sample_dir"/*_unaligned.fa; do

    [[ -f "\$fa" ]] || continue

    unaligned_count=\$((unaligned_count + 1))

    filename=\$(basename "\$fa")
    base=\${filename%_unaligned.fa}

    # Gene comes from the FASTA filename.
    gene=\${base%%_*}

    # Retrieve the matching canonical CDS 15-mer ONCE for this gene/file.
    cds_seq=\$(get_cds "\$gene")
    cds_seq=\$(printf '%s' "\$cds_seq" | tr '[:lower:]' '[:upper:]')
    canonical_stop_pos=\$(get_canonical_stop_pos "\$cds_seq")

    output_file="\$sample_out/\${base}_unaligned_categorized.tsv"

    echo
    echo "Processing UNALIGNED FASTA"
    echo "File: \$filename"
    echo "Gene: \$gene"
    echo "CDS 15-mer length: \${#cds_seq}"
    echo "Canonical CDS 15-mer STOP pos: \$canonical_stop_pos"
    echo "Output: \$output_file"

    echo -e "header\tatg1_pos\tatg0_pos\tdiff\tframe\tcanonical_cds_stop_pos\tatg1_frame_cds_stop_pos\tcategory" \
        > "\$output_file"

    awk '
        /^>/ {
            if (header != "")
                print header "\t" seq

            header=\$0
            seq=""
            next
        }
        {
            gsub(/[[:space:]\r]/, "", \$0)
            seq=seq \$0
        }
        END {
            if (header != "")
                print header "\t" seq
        }
    ' "\$fa" |
    while IFS=\$'\t' read -r header seq; do

        classify_sequence "\$header" "\$seq" "\$gene" "\$cds_seq" "\$canonical_stop_pos" "\$output_file"

    done

done

echo
echo "Unaligned FASTA files processed: \$unaligned_count"

# Clean up

rm -f "\$cds_lookup"

echo
echo "SLURM JOB FINISHED"
echo "Sample: \$sample"
echo "Date: \$(date)"

EOF


done

echo
echo "========================================"
echo "All samples submitted: $sample_count"
echo "========================================"