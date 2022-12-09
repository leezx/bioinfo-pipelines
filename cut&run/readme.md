
# 流程的执行步骤

- 第一步：single_chip_alignment_PE.sh 
- 第二步：single_chip_qc.sh
- 第三步：single_chip_misc_PE.sh
- 第四步：makemapper.sh
- 第五步：chip_group_qc.sh


## 第一步：single_chip_alignment_PE.sh 

```
trim_galore --paired --fastqc --stringency 3 --length 25 ${file} ${file/_1.fastq/_2.fastq}

bowtie2 --no-unal -p 5 -t -x /mnt/storage/home/da528/ref/hg19/hg19 -1 $file -2 $file2 -S ${file%.fastq}.sam 2>> ${file%.fastq}.$rand.log

samtools view -b -S ${file%.fastq}.unique.sam > ${file%.fastq}.unique.bam 2>> ${file%.fastq}.$rand.log

samtools sort -@ 3 ${file%.fastq}.unique.bam -o ${file%.fastq}.unique.sorted.bam 2>> ${file%.fastq}.$rand.log

samtools sort -@ 3 -n ${file%.fastq}.unique.sorted.bam | samtools view -bf 0x2 - | bedtools bamtobed -bedpe -i - 2> /dev/null | awk -v OFS="\t" '$1!="chrM" && $4!="chrM"{print $1,$2,$4,$6,$9,$10}' | sort | uniq -c | awk -v OFS="\t" '$1==1{m1++} $1==2{m2++} {md++} {mt+=$1} END {if(mt==0){mdt=-1} else {mdt=md/mt} if(md==0) {m1d=-1} else {m1d=m1/md} if(m2==0) {m12=-1} else {m12=m1/m2} printf "Total Uniquely Mapped Reads:\t%d\nTotal Distinct Genomic Locations:\t%d\nTotal Non-duplicated Reads:\t%d\nTotal 2-plicated reads:\t%d\nNRF:\t%f\nPBC1:\t%f\nPBC2:\t%f",mt,md,m1,m2,mdt,m1d,m12}' > qc/PBC_NRF/${file/.fastq/_PBC_NRF.txt}

samtools rmdup ${file%.fastq}.unique.sorted.bam ${file%.fastq}.bowtie2.bam 2>> ${file%.fastq}.$rand.log

samtools index $bam
```

## 第二步：single_chip_qc.sh

```
bam_stat.py -i $file > qc/base-stats/$base.base-stats.txt

read_NVC.py -i $file -o qc/NVC/$base

read_GC.py -i $file -o qc/GC/$base

Rscript run_spp.R -rf -odir=qc/cross-correlation -p=2 -savp -out=qc/cross-correlation/$base.log -c=$file

```

## 第三步：single_chip_misc_PE.sh

```
fastqc --noextract -j java -o qc/fastqc $file

macs callpeak -t $file -f BAMPE -n qc/macs/$base.macs14 -g $shortgenome

macs2 callpeak -t $file -f BAMPE -n qc/macs/$base.macs2 -q 0.01 -g $shortgenome

bamCoverage --binSize 25 -p 4 --normalizeUsing RPKM --extendReads $fragLength -b $file -of bigwig -o qc/bw/$base.bw

bamCoverage --binSize 10 --smoothLength 30 -p 4 --normalizeUsing RPGC --effectiveGenomeSize 2451960000 --extendReads $fragLength -b $file -of bigwig -o qc/bw/$base.10x30.bw
```

## 第四步：makemapper.sh

```
#
```

## 第五步：chip_group_qc.sh

```
#
```

