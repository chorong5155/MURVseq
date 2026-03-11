#!/bin/bash
account="p31754"
queue="long"
genomeID="hg38+VacV"
refDir="/projects/b1154/referenceGenomes"
refSeq="${refDir}/${genomeID}/${genomeID}.fna"
topDir=$(pwd)
outDir="$topDir/debug"
fastqDir="$topDir/fastq"

mkdir -p "$outDir"

ppnAlign=24


for umiR1 in "$fastqDir"/*_umi_1P.fastq.gz; do
    name=$(basename "${umiR1%_umi_1P.fastq.gz}" )
    umiR2="${umiR1/_1P/_2P}"
    outputdir_star="${topDir}/output/${name}"

    mkdir -p "$outputdir_star"

    sbatch <<EOF
#!/bin/bash -l
#SBATCH -A $account
#SBATCH -p $queue
#SBATCH -o $outDir/STAR${name}.out
#SBATCH -e $outDir/STAR${name}.err
#SBATCH -t 4:00:00
#SBATCH -c ${ppnAlign}
#SBATCH --mem=100g
#SBATCH -J STAR_${name}


echo "Running STAR alignment for sample: $name"
date

conda activate tt-seq

STAR --runThreadN ${ppnAlign} \
--runMode alignReads \
--genomeDir ${refDir}/${genomeID} \
--readFilesIn $umiR1 $umiR2 \
--readFilesCommand gunzip -c \
--outFilterType BySJout \
--outFilterMultimapNmax 20 \
--alignSJoverhangMin 1 \
--alignSJDBoverhangMin 1 \
--outFilterMismatchNmax 999 \
--outFilterMismatchNoverReadLmax 0.33 \
--alignIntronMin 20 \
--alignIntronMax 1000000 \
--alignMatesGapMax 1000000 \
--quantMode TranscriptomeSAM GeneCounts \
--outReadsUnmapped Fastx \
--outSAMtype BAM SortedByCoordinate \
--alignEndsType EndToEnd \
--limitBAMsortRAM 150000000000 \
--outSAMmultNmax 1 \
--outFileNamePrefix $outputdir_star/${name}_

# Post-processing
samtools index ${outputdir_star}/${name}_Aligned.sortedByCoord.out.bam

samtools view -@ ${ppnAlign} -b -f 2 ${outputdir_star}/${name}_Aligned.sortedByCoord.out.bam > ${outputdir_star}/${name}_Aligned_properlyPaired.sortedByCoord.out.bam
samtools index ${outputdir_star}/${name}_Aligned_properlyPaired.sortedByCoord.out.bam

# Deduplication
STAR --runThreadN ${ppnAlign} \
--runMode inputAlignmentsFromBAM \
--inputBAMfile ${outputdir_star}/${name}_Aligned_properlyPaired.sortedByCoord.out.bam \
--bamRemoveDuplicatesType UniqueIdentical \
--outSAMtype BAM SortedByCoordinate \
--limitBAMsortRAM 100000000000 \
--outFileNamePrefix ${outputdir_star}/${name}_markdups_

samtools view -@ ${ppnAlign} -b -F 1024 ${outputdir_star}/${name}_markdups_Processed.out.bam > ${outputdir_star}/${name}_nodups_Processed.out.bam
samtools index ${outputdir_star}/${name}_nodups_Processed.out.bam

# Generate signal tracks
STAR --runThreadN ${ppnAlign} \
--runMode inputAlignmentsFromBAM \
--inputBAMfile ${outputdir_star}/${name}_nodups_Processed.out.bam \
--outWigType bedGraph read1_5p \
--outWigStrand Unstranded \
--outWigNorm RPM \
--limitBAMsortRAM 100000000000 \
--outFileNamePrefix ${outputdir_star}/${name}_nodups_tracks_

date
echo "Finished STAR pipeline for: $name"
EOF

done
