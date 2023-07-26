#!/bin/bash

#This script is designed to take Illumina paired end read fq.gz files and process them using trimmomatic, STAR and HTSeq Union
#Best practise Illumina RNAseq processing steps found in Luis A. Corchete et al 2020
#other details about the RNAseq
#library prep was unstranded (htseq-count -s no )
#illumina device was V3 hiseq or newer therefore using Truseq3-PE-2 adapters (trimmomatic)
#reference genome = Homo_sapiens.GRCh38.dna.primary_assembly.fa (star)
#genome annotation = Homo_sapiens.GRCh38.109.gtf (star)

# Display Help
help()
{
  echo "This script is designed to take Illumina paired end read fq.gz files and process them using trimmomatic, RUM, HTSeq Union and TMM"
  echo "Raw data files (*.fq.gz) from Illumina runs should be placed into a custom study directory"
  echo
  echo "Syntax: scriptTemplate [-h|-r]"
  echo "options:"
  echo "h     Print this help"
  echo "r     Raw data directory"
  echo "e     conda environment defaults to illumina23"
  echo
}

#set options
conda_environment="illumina23"

#look for flags
while getopts ":h:r:e:t:" opt; do
  case $opt in
    h) help
    ;;
    r) rawdata_dir=$OPTARG
    r_option_present=true
    ;;
    e) conda_environment=$OPTARG
    ;;
    t) threads=$OPTARG
    ;;
    ?) echo "Invalid option: -$OPTARG"
            exit 1
    ;;
  esac
done

if [ "$1" == "-h" ]; then
  help
  exit 0
fi

echo -e "\n### 1. Raw data ###\n"
#Is rawdata_dir provided
if [ ! -d "$rawdata_dir" ]; then
  echo -e "\nError: valid raw data directory not specified - please supply using -r \n"
  exit 1
else
  echo -e "\nRaw data directory : $rawdata_dir is valid \n"
fi

# Are raw data .fq.qz files in rawdata_dir
if [ -n "$(find "$rawdata_dir" -name "*.fq.gz" -print -quit)" ]; then
  echo -e "\nRaw files '.fq.gz' found in the '$rawdata_dir' \n"
else
  echo -e "\nError: no raw files '.fq.gz' found in the '$rawdata_dir' \n"
  exit 1
fi

echo -e "\n### 2. Activate Conda env ### \n"
# Activate conda_environment
eval "$(conda shell.bash hook)"
conda activate $conda_environment
echo -e "\nActivated conda env: '$conda_environment' \n"
# check conda environment
# conda info --env

# fastQC
echo -e "\n### 3. FastQC ### \n"

# make a directory for fastqc
fastqc_dir="$(dirname "$rawdata_dir")/03.fastqc_reports"
if [ ! -d "$fastqc_dir" ]; then
    mkdir -p "$fastqc_dir"
    echo -e "\nCreating trim_output directory '$fastqc_dir' \n"
  else
    echo -e "\nDirectory: '$fastqc_dir' already exists \n"
fi

# Check if fastqc has been performed (maxdepth 1 = only first level of directorys, -print prints the file name -quit only finds the first file and then quits the search)
if [ -n "$(find "$fastqc_dir" -maxdepth 1 -name "*.html" -print -quit)" ]; then
  echo -e "\nSkipping fastQC as fastqc_reports directory '$fastqc_dir' contains fastqc html files \n"
else
  echo -e "\nStarting fastQC \n"
  read1_pattern="*_1.fq.gz"
  sub1="_1\.fq\.gz"
  sub2="_2.fq.gz"
  # find paired-end FASTQ files and run FastQC on them
  find "$rawdata_dir" -type f -name "$read1_pattern" -exec sh -c '
      read1="$1"
      read2=$(echo "$read1" | sed "s/'$sub1'/'$sub2'/")

      echo "Read 1: $read1"
      echo "Read 2: $read2"
      #echo "output: '$fastqc_dir'"
      fastqc "$read1" "$read2" -o '$fastqc_dir'
  ' sh {} \;
fi

# multiqc of fastqc RawData
if [ -n "$(find "$fastqc_dir" -maxdepth 1 -name "multiqc_report*.html" -print -quit)" ]; then
  echo -e "\nSkipping mulitQC as fastqc_reports directory '$fastqc_dir' contains a multiqc html file \n"
else
  echo -e "\nStarting multiQC \n"
  multiqc "$fastqc_dir"/*fastqc.zip -o "$fastqc_dir"
fi

echo -e "\n### 4. Trim ### \n"

# Make trimmomatic output folder
trim_dir="$(dirname "$rawdata_dir")/04.trim_output"
if [ ! -d "$trim_dir" ]; then
    mkdir -p "$trim_dir"
    echo -e "\nCreating trim_output directory '$trim_dir' /n"
  else
    echo -e "\nDirectory: '$trim_dir' already exists \n"
fi

# Check if files have already been trimmed
if [ -n "$(find "$trim_dir" -maxdepth 1 -name "*.fastq.gz" -print -quit)" ]; then
echo -e "\nSkipping trim as trim directory '$trim_dir' contains .fastqc.gz files \n"
else
  echo -e "\nstarting trim \n"
  read1_pattern="*_1.fq.gz"
  sub1="_1\.fq\.gz"
  sub2="_2.fq.gz"
  sub1_trim="_1_trimmed.fastq.gz"
  sub2_trim="_2_trimmed.fastq.gz"
  sub1_unpaired="_1_unpaired.fastq.gz"
  sub2_unpaired="_2_unpaired.fastq.gz"

  #find fq.gz paired reads and trim together
  find "$rawdata_dir" -type f -name "$read1_pattern" -exec sh -c '
      read1="$1"
      #passed $trim_dir as an argument to the subshell in find command at the bottom
      trim_dir="$2"
      #substitute 1 for 2 in the read name
      read2=$(echo "$read1" | sed "s/'$sub1'/'$sub2'/")
      #substitute read1 to make file names for trimmed and unpaired files and only take the basename.
      trim1_basename=$(basename $(echo "$read1" | sed "s/'$sub1'/'$sub1_trim'/"))
      trim2_basename=$(basename $(echo "$read1" | sed "s/'$sub1'/'$sub2_trim'/"))
      unpaired1_basename=$(basename $(echo "$read1" | sed "s/'$sub1'/'$sub1_unpaired'/"))
      unpaired2_basename=$(basename $(echo "$read1" | sed "s/'$sub1'/'$sub2_unpaired'/"))
      #add the trim directory to each basename so they are saved in trim_dir
      trim1="$2/$trim1_basename"
      trim2="$2/$trim2_basename"
      unpaired1="$2/$unpaired1_basename"
      unpaired2="$2/$unpaired2_basename"

      #echo "Read 1: $read1"
      #echo "Read 2: $read2"
      #echo "Trim 1: $trim1"
      #echo "Trim 2: $trim2"
      #echo "Unpaired 1: $unpaired1"
      #echo "Unpaired 2: $unpaired2"

      threads="$3"

      #Run trim
      trimmomatic PE -threads "$threads" "$read1" "$read2" "$trim1" "$unpaired1" "$trim2" "$unpaired2" ILLUMINACLIP:/home/minion-analysis/anaconda3/envs/illumina23/share/trimmomatic-0.39-2/adapters/TruSeq3-PE-2.fa:2:30:10:7:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:51
  ' sh {} "$trim_dir" "$threads" \;
fi

echo -e "\n### 5. Trim FastQC ### \n"

# make a directory for fastqc for trimmed reads
trim_fastqc_dir="$(dirname "$rawdata_dir")/05.trim_fastqc_reports"
if [ ! -d "$trim_fastqc_dir" ]; then
    mkdir -p "$trim_fastqc_dir"
    echo -e "\nCreating trim_output directory '$trim_fastqc_dir' \n"
  else
    echo -e "\nDirectory: '$trim_fastqc_dir' already exists \n"
fi

# Check if fastqc has been performed (maxdepth 1 = only first level of directorys, -print prints the file name -quit only finds the first file and then quits the search)
if [ -n "$(find "$trim_fastqc_dir" -maxdepth 1 -name "*.html" -print -quit)" ]; then
  echo -e "\nSkipping trimmed read fastQC as trim_fastqc_reports directory '$trim_fastqc_dir' contains fastqc html files \n"
else
  echo -e "\nStarting trimmed read fastQC \n"
  trim1_pattern="*_1_trimmed.fastq.gz"
  sub1="_1_trimmed\.fastq\.gz"
  sub2="_2_trimmed.fastq.gz"
  # find paired-end FASTQ files and run FastQC on them
  find "$trim_dir" -type f -name "$trim1_pattern" -exec sh -c '
      trim1="$1"
      trim2=$(echo "$trim1" | sed "s/'$sub1'/'$sub2'/")

      echo "Trimmed read 1: $trim1"
      echo "Trimmed read 2: $trim2"
      #echo "output: '$trim_fastqc_dir'"
      fastqc "$trim1" "$trim2" -o '$trim_fastqc_dir'
  ' sh {} \;
fi

# multiqc of trimmed fastqc RawData
if [ -n "$(find "$trim_fastqc_dir" -maxdepth 1 -name "multiqc_report*.html" -print -quit)" ]; then
  echo -e "\nSkipping trimmed read mulitQC as trim_fastqc_reports directory '$trim_fastqc_dir' contains a multiqc html file \n"
else
  echo -e "\nStarting multiQC \n"
  multiqc "$trim_fastqc_dir"/*fastqc.zip -o "$trim_fastqc_dir"
fi

echo -e "\n### 6. Align with STAR ### \n"
# STAR
# make a directory
align_dir="$(dirname "$rawdata_dir")/06.align_output"
if [ ! -d "$align_dir" ]; then
    mkdir -p "$align_dir"
    echo -e "\nCreating directory '$align_dir' \n"
  else
    echo -e "\nDirectory: '$align_dir' already exists \n"
fi

#check if alignment has been done then perform STAR alignment using the GRCh38 genome index in star_genome_indices. 3 bam outputs: aligned = genome, aligned.sortedbycoord = genome but sorted and aligned.toTranscriptome = genome alignment projected onto transcriptome (Not actually using this aligned to transcriptome file).
#Also need to index the aligned bam files and save into the save folder with the same name (but ending in .bai) for htseq count in the next step
if [ -n "$(find "$align_dir" -maxdepth 1 -name "*.bam" -print -quit)" ]; then
  echo -e "\nSkipping alignment as directory: '$align_dir' contains bam files \n"
else
  echo -e "\nStarting alignment \n"
  trim1_pattern="*_1_trimmed.fastq.gz"
  sub1="_1_trimmed\.fastq\.gz"
  sub2="_2_trimmed.fastq.gz"
  transcriptome_align_suffix="Aligned.toTranscriptome.out.bam"
  transcriptome_align_sorted_suffix="Aligned.toTranscriptome.sortedByCoord.out.bam"
  genome_align_suffix="Aligned.sortedByCoord.out.bam"
  # find paired-end FASTQ files and run FastQC on them
  find "$trim_dir" -type f -name "$trim1_pattern" -exec sh -c '
      trim1="$1"
      trim2=$(echo "$trim1" | sed "s/'$sub1'/'$sub2'/")
      trim1_basename=$(basename "$trim1")
      align_prefix="$2/${trim1_basename%'$sub1'}_"
      transcriptome_align="$align_prefix"'$transcriptome_align_suffix'
      transcriptome_align_sorted="$align_prefix"'$transcriptome_align_sorted_suffix'
      genome_align="$align_prefix"'$genome_align_suffix'

      echo "Trimmed read 1: $trim1"
      echo "Trimmed read 2: $trim2"
      echo "output: $align_prefix"
      threads="$3"

      STAR --runThreadN "$threads" --genomeDir /mnt/ssd2/star_genome_indices --readFilesIn "$trim1" "$trim2" --readFilesCommand gunzip -c --quantMode TranscriptomeSAM GeneCounts --outFileNamePrefix "$align_prefix" --outSAMtype BAM Unsorted SortedByCoordinate

      #echo "Sorting transcriptome alignment"
      #Not using transcriptome alignment file so no need to include the lengthy process of sorting it
      #echo "transcriptome aligned name: $transcriptome_align"
      #echo "transcriptome aligned sorted name: $transcriptome_align_sorted"
      #samtools sort "$transcriptome_align" -o "$transcriptome_align_sorted"

      echo "Indexing bam file: $genome_align"
      samtools index "$genome_align"
  ' sh {} "$align_dir" "$threads" \;
fi

echo -e "\n### 7. Aligned QC - flagstat ### \n"
# make a directory
align_qc_dir="$(dirname "$rawdata_dir")/07.align_QC"
if [ ! -d "$align_qc_dir" ]; then
    mkdir -p "$align_qc_dir"
    echo -e "\nCreating directory '$align_qc_dir' \n"
  else
    echo -e "\nDirectory: '$align_qc_dir' already exists \n"
fi

#check if flagstat has been performed and then perform flagstat on aligned BAM files
if [ -n "$(find "$align_qc_dir" -maxdepth 1 -name "*.flagstat" -print -quit)" ]; then
  echo -e "\nSkipping alignment QC as directory: '$align_qc_dir' contains flagstat files \n"
else
  echo -e "\nStarting alignment QC \n"
  align_pattern="*.out.bam"
  # find paired-end FASTQ files and run FastQC on them
  find "$align_dir" -type f -name "$align_pattern" -exec sh -c '
      align1="$1"
      align1_basename=$(basename "$align1")
      align1_output="$2/$align1_basename.flagstat"


      echo "Align file: $align1"
      echo "Output file: $align1_output"
      echo "Running samtools flagstat"
      samtools flagstat "$align1" > "$align1_output"

  ' sh {} "$align_qc_dir" \;
fi

echo -e "\n### 8. Count ### \n"
# HTSeq Union
count_suffix="gene_count.txt"
count_dir="$(dirname "$rawdata_dir")/08.count_output"
if [ ! -d "$count_dir" ]; then
    mkdir -p "$count_dir"
    echo -e "\nCreating directory '$count_dir' \n"
  else
    echo -e "\nDirectory: '$count_dir' already exists \n"
fi

if [ -n "$(find "$count_dir" -maxdepth 1 -name "*.txt" -print -quit)" ]; then
  echo -e "\nSkipping alignment as directory: '$count_dir' contains txt files \n"
else
  echo -e "\nStarting gene count \n"
  genome_align_pattern="*Aligned.sortedByCoord.out.bam"
  genome_align_suffix="Aligned.sortedByCoord.out.bam"
  gtf="/mnt/ssd2/reference/Homo_sapiens.GRCh38.109.gtf"
  # find paired-end FASTQ files and run FastQC on them
  find "$align_dir" -type f -name "$genome_align_pattern" -exec sh -c '
      align1="$1"
      align1_basename=$(basename "$align1")
      count_file="${align1_basename%'$genome_align_suffix'}'$count_suffix'"
      count="$2/$count_file"
      gtf='$gtf'

      echo "Alignment file: $align1"
      echo "Count output: $count"
      echo "GTF path: $gtf"
      echo "Running htseq-count"
      htseq-count -m union -s no -r pos -a 10 -i gene_id -f bam "$align1" "$gtf" > "$count"
  ' sh {} "$count_dir" \;
fi

echo -e "\n### 9. Count QC ### \n"
#create bash script to summarise the counting results
countqc_dir="$(dirname "$rawdata_dir")/09.count_QC"
if [ ! -d "$countqc_dir" ]; then
    mkdir -p "$countqc_dir"
    echo -e "\nCreating directory '$countqc_dir' \n"
  else
    echo -e "\nDirectory: '$countqc_dir' already exists \n"
fi

if [ -n "$(find "$countqc_dir" -maxdepth 1 -name "*.txt" -print -quit)" ]; then
  echo -e "\nSkipping alignment as directory: '$countqc_dir' contains txt files \n"
else
  echo -e "\nStarting gene count QC \n"
  # MulitQC for the Htseq count results
  multiqc "$count_dir"/*gene_count.txt -o "$countqc_dir"
  count_pattern="*${count_suffix}"

  # split the count file into head (genes with reads) and tail (summary of not aligned reads) and then to summarise the gene data with how many have reads and how many don't
  find "$count_dir" -type f -name "$count_pattern" -exec sh -c '
      count1="$1"
      echo "count file: $count1"
      count1_basename=$(basename "$count1")
      count_suffix="$2"
      #remove suffix
      countqc_file="${count1_basename%$count_suffix}"
      #echo "count QC file: $countqc_file"
      #add countqc_dir path
      count1qc="$3/$countqc_file"
      #new output files
      head_suffix="gene_count_head.txt"
      tail_suffix="gene_count_tail.txt"
      summary_suffix="gene_count_summary.txt"
      count1qc_head="$count1qc$head_suffix"
      #echo "head file: $count1qc_head"
      count1qc_tail="$count1qc$tail_suffix"
      #echo "tail file: $count1qc_tail"
      count1qc_summary="$count1qc$summary_suffix"
      #echo "summary file: $count1qc_summary"

      #Save the head (gene count data) and tail (unaligned reads summary)
      head -n -5 $count1 > $count1qc_head
      tail -n 5  $count1 > $count1qc_tail

      # Count the number of genes with expression greater than 0
          genes_with_expression_gt_0=$(awk '\''$2 > 0 { count++ } END { print count }'\'' "$count1qc_head")
          # count the number of genes with expression greater than 10
          genes_with_expression_gt_10=$(awk '\''$2 >10 {count++} END {print count}'\'' "$count1qc_head")
          # Count the number of genes with expression equal to 0
          genes_with_expression_eq_0=$(awk '\''$2 == 0 { count++ } END { print count }'\'' "$count1qc_head")

          # Create a summary file
          echo "File: $count1" > "$count1qc_summary"
          echo "__genes_with_reads_>10 $genes_with_expression_gt_10" >> "$count1qc_summary"
          echo "__genes_with_reads_>0 $genes_with_expression_gt_0" >> "$count1qc_summary"
          echo "__genes_no_reads $genes_with_expression_eq_0" >> "$count1qc_summary"
          cat "$count1qc_tail" >> "$count1qc_summary"

  ' sh {} "$count_suffix" "$countqc_dir" \;
fi

#To do
#run this and then take the results to Rstudio to do limma
#limma should be used for DE but has to be done with R script
#normalisation in limma
# glm to add the RIN data
#save to github
echo -e "\n### . ### \n"

echo -e "\n### Fin ### \n"

conda deactivate
