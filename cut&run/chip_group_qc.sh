#! /bin/bash
cd $1
genome=$2
name=$3

theBams=$(echo $@ | sed 's/ /\n/g' | grep 'bowtie2.bam' | tr '\n' ' ')

#echo makemapper -g $genome -o $name -k 10 -i bams/*.bam $theBams -q -d $1
#run elsewhere
#sh /ifs/labs/shivdasani/sv681/ref/bin/makemapper.sh -g $genome -o $name -c 10 -i $theBams -q -d $1

#mkdir -p fakemaker_$3

#sh mapmaker_cluster.sh -g $genome -m 0.01 -q -n 1 -o try1 -s /ifs/labs/shivdasani/sv681/ref/$genome/$genome.promoter-1k2k.bed

theBases=$(echo $@ | sed 's/ /\n/g' | grep 'bowtie2.bam' | sed 's/bowtie2.bam$/bowtie2.base-stats.txt/g' | awk '{print "qc/base-stats/"$0}' | tr '\n' ' ')
thePBCs=$(echo $@ | sed 's/ /\n/g' | grep 'bowtie2.bam' | sed 's/.bowtie2.bam$/_PBC_NRF.txt/g' | awk '{print "qc/PBC_NRF/"$0}' | tr '\n' ' ')
theCCs=$(echo $@ | sed 's/ /\n/g' | grep 'bowtie2.bam' | sed 's/bowtie2.bam$/bowtie2.log/g' | awk '{print "qc/cross-correlation/"$0}' | tr '\n' ' ')
echo "running base_stats as group on" $theBams
awk 'NF==2&&ARGIND==2{print $1}' FS=":" $theBases > qc/base-stats/$3.base-stats

#pull together bam-stats
for var in $theBases
do 
    printf "\t"$var >> qc/base-stats/$3.header.txt
    awk 'NF==2{print $2}' FS=":" $var | sed 's/ //g' | paste qc/base-stats/$3.base-stats - > qc/base-stats/$3.base-stats2
    mv qc/base-stats/$3.base-stats2 qc/base-stats/$3.base-stats
done
echo "" >> qc/base-stats/$3.header.txt
sed 's/qc\/base-stats\///g' qc/base-stats/$3.header.txt | sed 's/bowtie2.base-stats.txt//g' | cat - qc/base-stats/$3.base-stats > qc/base-stats/$3.base-stats.txt
rm qc/base-stats/$3.header.txt
rm qc/base-stats/$3.base-stats
mv qc/base-stats/$3.base-stats.txt qc/

#Calculate frip score. NOTE: make sure peak file name is correct in relation to BAM file name.
echo -e "FRiP_score\nReads_in_Peaks\nTotal_records" > qc/frip_scores.txt
for var in $theBams 
do 
    numReads=$(samtools view -c $var)
    samtools view $var | awk '{print $3,$4,$4+length($10)}' OFS="\t" | sort -k 1,1 -k2,2n | bedtools coverage -counts -sorted -a qc/macs/${var/.bam/.macs*_peaks.narrowPeak} -b - | awk 'BEGIN {sum=0} {sum=sum+$NF} END {print sum/'$numReads'"\n"sum"\n"'$numReads'}' | paste qc/frip_scores.txt - > qc/temp
    mv qc/temp qc/frip_scores.txt
done

#pull together the PBC/NRF score stuff
awk 'ARGIND==2{print $1}' FS=":" $thePBCs > qc/PBCs.txt
for var in $thePBCs
do 
    awk '{print $2}' FS=":\t" $var | paste qc/PBCs.txt - > qc/temp
    mv qc/temp qc/PBCs.txt
done

#pull together the cross-correlation stuff
awk 'FNR==1&&ARGIND==2{print $2"\n"$3"\n"$4"\n"$5"\n"$6"\n"$7"\n"$8"\n"$9"\n"$10"\n"$11}' FS="\t" $theCCs > qc/CCs.txt
for var in $theCCs
do 
    awk 'NR==2{print $2"\n"$3"\n"$4"\n"$5"\n"$6"\n"$7"\n"$8"\n"$9"\n"$10"\n"$11}' $var | paste qc/CCs.txt - > qc/temp
    mv qc/temp qc/CCs.txt
done

cat qc/$3.base-stats.txt qc/frip_scores.txt qc/PBCs.txt qc/CCs.txt benchmarks.txt > qc/summary.tsv

#theBases=$(echo $@ | sed 's/ /\n/g' | grep 'fastq' | sed 's/fastq.gz$/STAR-encode.read_dist.txt/g' | sed 's/fastq$/STAR-encode.read_dist.txt/g' | awk '{print "rnaqc/read_dist/"$0}' | tr '\n' ' ')
#echo "running read_dist as group on" $theBams
#awk 'NF>3&&ARGIND==1&&NR!=5{print $1}' $theBases > rnaqc/read_dist/$name.read_dist; for var in $theBases; do printf "\t"$var >> rnaqc/read_dist/$name.header.txt; awk 'NF>3&&NR!=5{print $4}' $var | paste rnaqc/read_dist/$name.read_dist - > rnaqc/read_dist/$name.read_dist2; mv rnaqc/read_dist/$name.read_dist2 rnaqc/read_dist/$name.read_dist; done; echo "" >> rnaqc/read_dist/$name.header.txt; sed 's/rnaqc\/read_dist\///g' rnaqc/read_dist/$name.header.txt | sed 's/STAR-encode.read_dist.txt//g' | cat - rnaqc/read_dist/$name.read_dist > rnaqc/read_dist/$name.read_dist.txt; rm rnaqc/read_dist/$name.header.txt; rm rnaqc/read_dist/$name.read_dist

#theBases=$(echo $@ | sed 's/ /\n/g' | grep 'fastq' | sed 's/fastq.gz$/STAR-encode.counts/g' | sed 's/fastq$/STAR-encode.counts/g' | awk '{print "rnaqc/HTseq/"$0}' | tr '\n' ' ')
#echo "running htseq as group on" $theBams
#awk 'substr($1,1,2)!="__"&&NR!=1&&ARGIND==2{print $1}' FS=":" $theBases > rnaqc/HTseq/$3.counts; for var in $theBases; do printf "\t"$var >> rnaqc/HTseq/$3.header.txt; awk 'substr($1,1,2)!="__"&&NR!=1{print $4}' FS=":" $var | sed 's/ //g' | paste rnaqc/HTseq/$3.counts - > rnaqc/HTseq/$3.counts2; mv rnaqc/HTseq/$3.counts2 rnaqc/HTseq/$3.counts; done; echo "" >> rnaqc/HTseq/$3.header.txt; sed 's/rnaqc\/HTseq\///g' rnaqc/HTseq/$3.header.txt | sed 's/STAR-encode.counts//g' | cat - rnaqc/HTseq/$3.counts > rnaqc/HTseq/$3.counts.tsv; rm rnaqc/HTseq/$3.header.txt; rm rnaqc/HTseq/$3.counts

#echo "creating deseq2 index on" $theBams
#printf "sample\tlocation\tall\ttype\n" > rnaqc/HTseq/$3.deseq2-index.tsv
#echo $@ | sed 's/ /\n/g' | grep 'fastq' | sed 's/fastq.gz$/STAR-encode.counts/g' | sed 's/fastq$/STAR-encode.counts/g' | awk '{printf "\t"$0"\tall\n"}' >> rnaqc/HTseq/$3.deseq2-index.tsv

#theBams=$(echo $@ | sed 's/ /\n/g' | grep 'fastq' | sed 's/fastq.gz$/STAR-encode.bam/g' | sed 's/fastq$/STAR-encode.bam/g' | tr '\n' ' ')
#mkdir -p rnaqc/genebody_cov
#echo "running genebody_coverage on" $theBams
#/ifs/labs/shivdasani/sv681/ref/python/geneBody_coverage.py -i $theBams -r /ifs/labs/shivdasani/sv681/ref/$genome/refSeq.bed -o rnaqc/genebody_cov/$name

