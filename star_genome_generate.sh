#!bin/bash

STAR --runThreadN 26 --runMode genomeGenerate --genomeDir /mnt/ssd2/star_genome_indices/ --genomeFastaFiles /mnt/ssd2/reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa --sjdbGTFfile /mnt/ssd2/reference/Homo_sapiens.GRCh38.109.gtf --sjdbOverhang 100
