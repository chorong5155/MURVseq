#!/bin/bash

emailaddress="${USER}@e.northwestern.edu"
topDir=$(pwd)
samplename=$(basename $topDir)
outputdir=${topDir}"/output"
genomeID="hg38+VacV"
scriptsDir="/projects/b1154/cpp4747/scripts"
trimmomaticPath="/home/cpp4747/software/anaconda3/share/trimmomatic-0.39-2"
fastq_screenDir="/home/cpp4747/software/anaconda3/share/fastq-screen-0.13.0-0/fastq_screen"
refDir="/projects/b1154/referenceGenomes"
refSeq="${refDir}/${genomeID}/${genomeID}.fna"
account="b1042"
queue="genomics"
queue_time="12:00:00"
ppnAlign=24
outDir="$topDir/debug"
read1str="_R1" 
read2str="_R2" 
fastqdir=${topDir}"/fastq/*_R*.fastq*"
outputdir=${topDir}"/output"
tmpdir=${topDir}"/tmp"

if [ ! -e "${refDir}/${genomeID}/SA" ]; then
    echo "***! Reference sequence $refSeq does not appear to have been indexed. Please run STAR genomeGenerate on this file before running."
    exit 100
fi

if [ -z "$genomePath" ]; then
    genomePath=$refDir/$genomeID/$genomeID.chrom.sizes
fi

if [ ! -d "$outDir" ]; then
    mkdir "$outDir"
    chmod 777 "$outDir"
fi

if [ ! -d "$tmpdir" ]; then
    mkdir "$tmpdir"
    chmod 777 "$tmpdir"
fi

if [ ! -d "$topDir/fastq" ]; then
    echo "Directory \"$topDir/fastq\" does not exist."
    echo "Create \"$topDir/fastq\" and put fastq files to be aligned there."
    exit 100
else 
    if stat -t ${fastqdir} >/dev/null 2>&1; then
        echo "(-: Looking for fastq files...fastq files exist"
    else
        echo "***! Failed to find any files matching ${fastqdir}"
        exit 1
    fi
fi

if ls $topDir/fastq/*trimmed* 1> /dev/null 2>&1; then
    echo "*** Trimmed fastq files exist, removing them..."
    trimmedFiles=$(ls $topDir/fastq/*trimmed*)
    for file in $trimmedFiles; do
        rm $file
    done
fi

if [[ -d "$outputdir" ]]; then
    echo "***! Move or remove directory \"$outputdir\" before proceeding."
    exit 100
else
    mkdir "$outputdir" || { echo "***! Unable to create ${outputdir}, check permissions." ; exit 100; }
fi

fastqR1files=${topDir}"/fastq/*_R1*.fastq*"

for read1file in ${fastqR1files}; do
    read2file="${read1file/_R1/_R2}"
    if [[ ! -f "$read2file" ]]; then
        echo "***! Cannot find read 2 for fastq $read1file. Please make sure read 2 exists before proceeding."
        exit 1
    else
        echo "Read 2 file for read 1 file $read1file is $read2file"
        name=$(basename "${read1file%$read1str*}" | cut -d'_' -f1-4)
        echo "Basename for these fastqs will be $name."

        
        trimR1files=($topDir/fastq/${name}_trimmed_1P.fastq.gz)
        trimR2files=($topDir/fastq/${name}_trimmed_2P.fastq.gz)

       
        sbatch <<EOF
#!/bin/bash -l
#SBATCH -A $account
#SBATCH -p $queue
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mem-per-cpu=4gb
#SBATCH -o $outDir/_qc_${name}-%j.out
#SBATCH -e $outDir/_qc_${name}-%j.err
#SBATCH -t 24:00:00
#SBATCH -J "qc_${name}"


module load fastqc/0.11.5
fastqc -t 4 $read1file --outdir=$topDir/qc_output
fastqc -t 4 $read2file --outdir=$topDir/qc_output

module load bwa
$fastq_screenDir --aligner bwa --subset 1000000 --force $read1file --outdir=$topDir/qc_output
$fastq_screenDir --aligner bwa --subset 1000000 --force $read2file --outdir=$topDir/qc_output

# downloaded human rRNA wget https://raw.githubusercontent.com/fallerlab/ARF/main/rRNAs/4V6X_human_rRNAs.fa -O human_rRNAs.fa

conda activate bbmap
/home/cpp4747/software/anaconda3/envs/bbmap/bbduk.sh -Xmx27g \
in1=$read1file \
in2=$read2file \
out1=$topDir/qc_output/${name}_unmatched_R1.fastq.gz \
out2=$topDir/qc_output/${name}_unmatched_R2.fastq.gz \
outm1=$topDir/qc_output/${name}_rRNA_R1.fastq.gz \
outm2=$topDir/qc_output/${name}_rRNA_R2.fastq.gz \
ref=/home/cpp4747/software/human_rRNAs.fa \
stats=$topDir/qc_output/${name}_rRNA_stats.txt
EOF
    fi
done


read1strP="_trimmed_1P"
read2strP="_trimmed_2P"

# Trim adaptor sequences from fastq

trimR1files=()
trimR2files=()
dependtrim="afterok"
for read1file in ${fastqR1files[@]}; do
sleep 0.5
name=$(basename ${read1file%$read1str*})
trimR1files+=($topDir/fastq/${name}_trimmed_1P.fastq.gz)
trimR2files+=($topDir/fastq/${name}_trimmed_2P.fastq.gz)

# trimmomaticVersion=$(trimmomatic -version)
trimParameters="ILLUMINACLIP:${trimmomaticPath}/adapters/enzyme_adaptors.fa:2:30:10:1:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50"

jid=`sbatch <<- TRIM | egrep -o -e "\b[0-9]+$"
#!/bin/bash -l
#SBATCH -A $account
#SBATCH -p $queue
#SBATCH -o $outDir/trim_${name}_${samplename}-%j.out
#SBATCH -e $outDir/trim_${name}_${samplename}-%j.err
#SBATCH -t 20:00:00
#SBATCH -n 1
#SBATCH -c ${ppnAlign}
#SBATCH --mem-per-cpu=4g
#SBATCH -J "trim_${name}_${samplename}"
sleep 1; chmod -R 777 $topDir/debug
date


echo "Trimming fastq files using trimmomatic from conda env (version 0.39) and trimming parameters $trimParameters"

trimmomatic PE -threads ${ppnAlign} -basein $read1file -baseout $topDir/fastq/${name}_trimmed.fastq.gz $trimParameters

date

TRIM`
dependtrim="$dependtrim:$jid"
done

trimR1list=$( IFS=$','; echo "${trimR1files[*]}" )
trimR2list=$( IFS=$','; echo "${trimR2files[*]}" )


#umi_tools

read1strU="_umi_1P"
read2strU="_umi_2P"

umiR1files=()
umiR2files=()
dependumi="afterok"
for read1file in ${trimR1files[@]}; do
sleep 0.5
files=($read1file)
R1=${files[$SLURM_ARRAY_TASK_ID]}
R2=${R1/_1P/_2P}

name=$(basename ${read1file%$read1strP*})

umiR1files+=($topDir/fastq/${name}_umi_1P.fastq.gz)
umiR2files+=($topDir/fastq/${name}_umi_2P.fastq.gz)

#UMI extract 
jid=`sbatch <<- UMI | egrep -o -e "\b[0-9]+$"
#!/bin/bash -l
#SBATCH -A $account
#SBATCH -p $queue
#SBATCH --array=6
#SBATCH -o $outDir/umi_${name}_${samplename}-%j.out
#SBATCH -e $outDir/umi_${name}_${samplename}-%j.err
#SBATCH -t 20:00:00
#SBATCH -n 1
#SBATCH --mem-per-cpu=4g
#SBATCH -d ${dependtrim}
#SBATCH -J "umi_${name}_${samplename}"
sleep 1; chmod -R 777 $topDir/debug
date

conda activate umi

umi_tools extract --extract-method=string --bc-pattern=NNNNNNNNNN --stdin="$R1" --read2-in="$R2" --stdout=${topDir}/fastq/${name}_umi_1P.fastq.gz --read2-out=${topDir}/fastq/${name}_umi_2P.fastq.gz


if [[ ! -f ${topDir}/fastq/${name}_umi_2P.fastq.gz ]]
then
	exit 100
else
	echo "Cho, =D ${topDir}/fastq/${name}_umi_2P.fastq.gz created successfully"
fi


date
UMI`
dependumi="$dependumi:$jid"

done



