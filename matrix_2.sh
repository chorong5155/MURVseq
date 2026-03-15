#!/bin/bash

# Paths
debug="/projects/b1042/WalshLab/cpp4747/19215FL-11/debug"
count_dir="/projects/b1042/WalshLab/cpp4747/19215FL-11/matrix"
out_file="${count_dir}/merged_5prime_matrix.txt"

# Create debug folder if not exists
mkdir -p "$debug"

# SLURM parameters
account="b1042"
queue="genomics"

# Submit the merge job
sbatch <<BATCH
#!/bin/bash -l
#SBATCH -A $account
#SBATCH -p $queue
#SBATCH -n 1
#SBATCH -c 8
#SBATCH --mem-per-cpu=8G
#SBATCH -t 24:00:00
#SBATCH -o ${debug}/merge_5prime-%j.out
#SBATCH -e ${debug}/merge_5prime-%j.err
#SBATCH -J merge_5prime_counts

cd "$count_dir"

echo "Starting merge at \$(date)"

awk '
BEGIN { OFS="\t" }

FNR==1 {
    sample = FILENAME
    sub(".*/","",sample)
    sub("_5prime_counts.txt","",sample)
    files[++n] = sample
}

{
    key = \$2 FS \$3
    data[key, sample] += \$1
    keys[key] = 1
}

END {
    printf "chr\tpos"
    for (i=1; i<=n; i++) printf "\t%s", files[i]
    print ""

    for (k in keys) {
        split(k,a,FS)
        printf "%s\t%s", a[1], a[2]
        for (i=1; i<=n; i++)
            printf "\t%d", data[k, files[i]]+0
        print ""
    }
}
' *_5prime_counts.txt > "$out_file"

echo "Merge completed at \$(date). Output saved to $out_file"
BATCH
