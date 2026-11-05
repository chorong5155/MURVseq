#!/bin/bash

topdir=$(pwd)
account="b1042"
queue="genomics"

BASE_DIR="/projects/b1042/WalshLab/cpp4747/19215FL-11/aligned_UTR"
gene_map="${topdir}/gene_15mer_map.tsv"

out_dir="${topdir}/uAUGtest"
DEBUG="${topdir}/debug"

mkdir -p "$out_dir" "$DEBUG"

echo "OUTER SCRIPT STARTED"
echo "topdir   = $topdir"
echo "BASE_DIR = $BASE_DIR"
echo "gene_map = $gene_map"
echo "out_dir  = $out_dir"
echo "DEBUG    = $DEBUG"

if [[ ! -d "$BASE_DIR" ]]; then
    echo "ERROR: BASE_DIR does not exist:"
    echo "$BASE_DIR"
    exit 1
fi

if [[ ! -f "$gene_map" ]]; then
    echo "ERROR: gene_map does not exist:"
    echo "$gene_map"
    exit 1
fi

echo "BASE_DIR and gene_map exist."

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

echo
echo "SLURM JOB STARTED"
echo "Sample: \$sample"
echo "Date: \$(date)"

if [[ ! -d "\$sample_dir" ]]; then
    echo "ERROR: sample directory does not exist:"
    echo "\$sample_dir"
    exit 1
fi

if [[ ! -f "$gene_map" ]]; then
    echo "ERROR: gene_map does not exist:"
    echo "$gene_map"
    exit 1
fi

mkdir -p "\$sample_out"


get_15mer() {
    local gene_phase="\$1"
    awk -v g="\$gene_phase" '
        BEGIN {
            FS="[[:space:]]+"
        }
        \$1 == g {
            print toupper(\$2)
            exit
        }
    ' "$gene_map"
}


classify_sequence() {

    local header="\$1"
    local seq="\$2"
    local gene_phase="\$3"
    local seq15="\$4"
    local output_file="\$5"

    seq=\$(printf '%s' "\$seq" | tr '[:lower:]' '[:upper:]')
    seq15=\$(printf '%s' "\$seq15" | tr '[:lower:]' '[:upper:]')
    seq_len=\${#seq}
    ATG0_POS=\$((seq_len + 1))
    
    
    ATG1_POS=\$(printf '%s' "\$seq" |
        grep -b -o "ATG" |
        head -n1 |
        cut -d: -f1)


    if [[ -z "\$ATG1_POS" ]]; then

        echo -e "\$header\tNA\t\$ATG0_POS\tNA\tNA\tNA\tNA\tCategory4_No_uAUG" \
            >> "\$output_file"

        return
    fi

    ATG1_POS=\$((ATG1_POS + 1))



    diff=\$((ATG0_POS - ATG1_POS))
    frame=\$((diff % 3))

    combined_seq="\${seq}\${seq15}"
    combined_len=\${#combined_seq}


   
    stop_found=0
    stop_pos="NA"
    stop_region="NA"
    stop_codon="NA"

    start0=\$((ATG1_POS - 1))

    for ((i=start0; i+2<combined_len; i+=3)); do

        codon=\${combined_seq:\$i:3}

        for stop in "\${STOP_CODONS[@]}"; do

            if [[ "\$codon" == "\$stop" ]]; then

                stop_found=1

                # Convert back to 1-based coordinate
                stop_pos=\$((i + 1))
                stop_codon="\$codon"

               
                if [[ "\$stop_pos" -lt "\$ATG0_POS" ]]; then

                    stop_region="UTR"

                else

                    stop_region="CDS15"

                fi

                break 2
            fi

        done
    done


  

    if [[ "\$stop_found" -eq 1 && "\$stop_region" == "UTR" ]]; then

        category="Category2_Upstream_Stopped"

        echo -e "\$header\t\$ATG1_POS\t\$ATG0_POS\t\$diff\t\$frame\t\$stop_pos\tNA\t\$category" \
            >> "\$output_file"

        return
    fi


   
    if [[ "\$frame" -eq 0 ]]; then
        if [[ "\$stop_found" -eq 1 && "\$stop_region" == "CDS15" ]]; then
            cds_stop_relative=\$((stop_pos - ATG0_POS + 1))
            category="Category1_Productive_Upstream"
        else
            cds_stop_relative="NA"
            category="Category5_No_CDS_Stop"

        fi


   
    else

        if [[ "\$stop_found" -eq 1 && "\$stop_region" == "CDS15" ]]; then
            cds_stop_relative=\$((stop_pos - ATG0_POS + 1))
            category="Category3_Overlapping_OutOfFrame"
        else

            cds_stop_relative="NA"
            category="Category5_No_CDS_Stop"

        fi

    fi



    echo -e "\$header\t\$ATG1_POS\t\$ATG0_POS\t\$diff\t\$frame\tNA\t\$cds_stop_relative\t\$category" \
        >> "\$output_file"
}


# PROCESS ALIGNED FASTA

aligned_count=0

for fa in "\$sample_dir"/*_highscore_aligned_seqs.fa; do
    [[ -f "\$fa" ]] || continue
    aligned_count=\$((aligned_count + 1))
    filename=\$(basename "\$fa")

    base=\${filename%_highscore_aligned_seqs.fa}
    gene_phase=\${base%_utr}
    seq15=\$(get_15mer "\$gene_phase")
    seq15=\$(printf '%s' "\$seq15" | tr '[:lower:]' '[:upper:]')

    output_file="\$sample_out/\${base}_highscore_aligned_seqs_categorized.tsv"

   
    echo "Processing ALIGNED FASTA"
    echo "File: \$filename"
    echo "Gene-phase: \$gene_phase"
    echo "15-mer: \$seq15"
    echo "15-mer length: \${#seq15}"
    echo "Output: \$output_file"

    if [[ -z "\$seq15" ]]; then

        echo "WARNING: No 15-mer found for:"
        echo "\$gene_phase"

    elif [[ \${#seq15} -ne 15 ]]; then

        echo "WARNING: 15-mer length is not 15 for \$gene_phase"
        echo "Sequence: \$seq15"
        echo "Length: \${#seq15}"

    fi


    echo -e "header\tatg1_pos\tatg0_pos\tdiff\tframe\tutr_stop_pos\tcds15_stop_pos\tcategory" \
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

        classify_sequence \
            "\$header" \
            "\$seq" \
            "\$gene_phase" \
            "\$seq15" \
            "\$output_file"

    done

done

echo
echo "Aligned FASTA files processed: \$aligned_count"



unaligned_count=0

for fa in "\$sample_dir"/*_unaligned.fa; do

    [[ -f "\$fa" ]] || continue

    unaligned_count=\$((unaligned_count + 1))

    filename=\$(basename "\$fa")

    base=\${filename%_unaligned.fa}


   
    gene_phase=\${base%_utr}

    seq15=\$(get_15mer "\$gene_phase")
    seq15=\$(printf '%s' "\$seq15" | tr '[:lower:]' '[:upper:]')

    output_file="\$sample_out/\${base}_unaligned_categorized.tsv"

    echo
    echo "Processing UNALIGNED FASTA"
    echo "File: \$filename"
    echo "Gene-phase: \$gene_phase"
    echo "15-mer: \$seq15"
    echo "15-mer length: \${#seq15}"
    echo "Output: \$output_file"


    if [[ -z "\$seq15" ]]; then

        echo "WARNING: No 15-mer found for:"
        echo "\$gene_phase"

    elif [[ \${#seq15} -ne 15 ]]; then

        echo "WARNING: 15-mer length is not 15 for \$gene_phase"
        echo "Sequence: \$seq15"
        echo "Length: \${#seq15}"

    fi


    echo -e "header\tatg1_pos\tatg0_pos\tdiff\tframe\tutr_stop_pos\tcds15_stop_pos\tcategory" \
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

        classify_sequence \
            "\$header" \
            "\$seq" \
            "\$gene_phase" \
            "\$seq15" \
            "\$output_file"

    done

done


echo
echo "Unaligned FASTA files processed: \$unaligned_count"

echo
echo "SLURM JOB FINISHED"
echo "Sample: \$sample"
echo "Date: \$(date)"

EOF

done

echo "All samples submitted: $sample_count"
