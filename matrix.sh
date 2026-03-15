#!/bin/bash
debug="/projects/b1042/WalshLab/cpp4747/19215FL-11/debug"
bam_dir="/projects/b1042/WalshLab/cpp4747/19215FL-11/output"
out_dir="/projects/b1042/WalshLab/cpp4747/19215FL-11/matrix"
mkdir -p "$out_dir" "$debug"

account="b1042"
queue="genomics"

# Step 1: Extract 5' positions for each BAM
for sampleDir in "$bam_dir"/*/; do
    sampleDirName=$(basename "$sampleDir")
    
    for bam in "$sampleDir"/*_Aligned_properlyPaired.sortedByCoord.out.bam; do
        bamSample=$(basename "$bam" .bam)
        echo "Submitting job for $bamSample"

        sbatch <<BATCH
#!/bin/bash -l
#SBATCH -A $account
#SBATCH -p $queue
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mem-per-cpu=4G
#SBATCH -o ${debug}/_matrix_${bamSample}-%j.out
#SBATCH -e ${debug}/_matrix_${bamSample}-%j.err
#SBATCH -t 24:00:00
#SBATCH -J "mat_${bamSample}"

module load samtools

samtools view "$bam" | awk '{
    if (and(\$2,16)) { print \$3 "\t" \$4+length(\$10)-1 }
    else { print \$3 "\t" \$4 }
}' | sort -k1,1 -k2,2n | uniq -c > "${out_dir}/${bamSample}_5prime_counts.txt"
BATCH

    done
done

