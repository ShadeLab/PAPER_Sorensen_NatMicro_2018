# Genome Assembly and Binning Workflow
## Software used in this Workflow
### [MEGAHIT](https://github.com/voutcn/megahit)
### [bbmap](https://jgi.doe.gov/data-and-tools/bbtools/)
### [samtools](https://github.com/samtools/samtools) (Version 1.3 or later)
### [MetaBAT](https://bitbucket.org/berkeleylab/metabat)
### [CheckM](http://ecogenomics.github.io/CheckM/)

## Assembly using MEGAHIT
```
megahit --12 combined.fastq.pe --k-list 27,37,47,57,67,77,87,97,107 -o Megahit_QC_Assembly/ -t $PBS_NUM_PPN
```
## Indexing Assembled contigs using bbmap
```
bbmap.sh ref=MA_Contigs.fa build=1 -Xmx215g
```
## Mapping reads to contigs using bbmap
```
bbmap.sh in=/mnt/research/ShadeLab/Sorensen/Cen01.anqdp.fastq.gz build=1 -Xmx215g out=Cen01_MA.sam
```
## Converting SAM files into BAM files and sorting them using SAMtools
```
samtools view -bS Cen01_MA.sam > Cen01_MA.bam
samtools sort -o Cen01_MA.bam -T Cen01_Sorted -@ 8 -m 8G Cen01_MA.bam
samtools index -b Sorted_BAM_files/*
```
## Genome Binning using METABAT
```
metabat -i ../Megahit_QC_Assembly/final.contigs.fa -v -a depth_20141007.txt -o METABAT_VerySpecific --saveTNF saved.tnf --saveDistance saved.dist -t 40 --veryspecific
```
## Quality Checking and Identification using CheckM
```
checkm lineage_wf -t 8 -x fa VerySpecific_Bins/ CheckM_VerySpecific_Output
checkm qa CheckM_VerySpecific_Output/lineage.ms CheckM_VerySpecific_Output/ > CheckM_VerySpecific_Results.txt
```
