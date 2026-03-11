#!/bin/bash

topdir=$(pwd)
account="b1042"
queue="genomicslong"
DEBUG="${topdir}/debug"
aligned_dir="${topdir}/aligned_UTR"
aligned_sort="${topdir}/aligned_sort"

mkdir -p "$DEBUG" "$aligned_sort"

for sample_dir in "$aligned_dir"/*; do
    [[ -d "$sample_dir" ]] || continue
    sample=$(basename "$sample_dir")

    echo "Submitting job for sample: $sample"
    
    mkdir -p "$aligned_sort/$sample"
    
    sbatch --export=ALL,sample_dir="$sample_dir",sample="$sample",aligned_sort="$aligned_sort" <<EOF
#!/bin/bash -l
#SBATCH -A $account
#SBATCH --partition=$queue
#SBATCH -J aligned_sort_${sample}
#SBATCH -o $DEBUG/aligned_sort_${sample}.%j.out
#SBATCH -e $DEBUG/aligned_sort_${sample}.%j.err
#SBATCH -t 10-00:00:00
#SBATCH -c 4
#SBATCH --mem-per-cpu=40G

MAXLEN=20

# Loop over all aligned FASTA for this sample
for f in \$sample_dir/*_highscore_aligned_seqs.fa; do
    base=\$(basename "\$f" .fa)
    short_out="\$aligned_sort/$sample/\${base}_short.fa"
    long_out="\$aligned_sort/$sample/\${base}_long.fa"

    awk -v maxlen="\$MAXLEN" '
        /^>/ {header=\$0; next}
        {
            seq=\$0
            if(length(seq) < maxlen) {
                print header
                print seq
            } else {
                print header > "'"\$long_out"'"
                print seq > "'"\$long_out"'"
            }
        }
    ' "\$f" > "\$short_out"
    echo "Created \$short_out and \$long_out"
done

# Loop over all unaligned FASTA for this sample
for f in \$sample_dir/*_unaligned.fa; do
    base=\$(basename "\$f" .fa)
    short_out="\$aligned_sort/$sample/\${base}_short.fa"
    long_out="\$aligned_sort/$sample/\${base}_long.fa"

    awk -v maxlen="\$MAXLEN" '
        /^>/ {header=\$0; next}
        {
            seq=\$0
            if(length(seq) < maxlen) {
                print header
                print seq
            } else {
                print header > "'"\$long_out"'"
                print seq > "'"\$long_out"'"
            }
        }
    ' "\$f" > "\$short_out"
    echo "Created \$short_out and \$long_out"
done
EOF

done
