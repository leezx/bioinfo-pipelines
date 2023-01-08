#! /bin/bash

cd $1

file=$2
genome=$3
rand=$4
trim=$5

if [ ${file: -3} == ".gz" ]
	then
		if [ ! -s $file ]
			then
				echo "file" $file "does not exist"
				exit 1
		fi
		if [ ! -s ${file%.gz} ]
			then
				echo "running unzip on" $file "and" ${file/_1/_2}
				gzip -d $file
				gzip -d ${file/_1.fastq/_2.fastq}
				file=${file%.gz}
			else
				echo "unzipped version of" $file "already exists"
		fi
fi 



if [ ${file: -6} == ".fastq" ]
	then
		if [ ! -s $file ]
			then
				echo "file" $file "does not exist"
				exit 1
		fi
		if [ $trim == "yes" ]
            then
                if [ ! -s ${file/.fastq/_val_1.fastq} ]
                    then
                        echo "running trim_galore on" $file
                        trim_galore --paired --fastqc --stringency 3 --length 25 ${file} ${file/_1.fastq/_2.fastq}
                        mkdir -p fastq
                        gzip $file
                        gzip ${file/_1.fastq/_2.fastq}
                        mv ${file}.gz fastq/
                        mv ${file/_1.fastq/_2.fastq.gz} fastq/
                        mkdir -p qc
                        mkdir -p qc/trimming_reports
                        mv ${file}_trimming_report.txt qc/trimming_reports
                        mv ${file/_1.fastq/_2.fastq}_trimming_report.txt qc/trimming_reports
                        mkdir -p qc/fastqc
                        mv ${file/.fastq/_val_1_fastqc}* qc/fastqc
                        mv ${file/_1.fastq/_2_val_2_fastqc}* qc/fastqc
                        mv ${file/.fastq/_val_1.fq} ${file/.fastq/_val_1.fastq}
                        mv ${file/_1.fastq/_2_val_2.fq} ${file/_1.fastq/_2_val_2.fastq}
                fi
                file=${file/.fastq/_val_1.fastq}
                file2=${file/_1_val_1/_2_val_2}
        fi
        if [ ! -s ${file/.fastq/.bowtie2.bam} ]
            then
                echo "running bowtie2 on" $file
                bowtie2 --no-unal -p 5 -t -x /mnt/storage/home/da528/ref/hg19/hg19 -1 $file -2 $file2 -S ${file%.fastq}.sam 2>> ${file%.fastq}.$rand.log
                
                #remove multimapping reads
                grep -v XS:i: ${file%.fastq}.sam > ${file%.fastq}.unique.sam 2>> ${file%.fastq}.$rand.log
                rm ${file%.fastq}.sam
                
                #convert to bam format
                samtools view -b -S ${file%.fastq}.unique.sam > ${file%.fastq}.unique.bam 2>> ${file%.fastq}.$rand.log
                rm ${file%.fastq}.unique.sam
                
                #sort the bam based on genomic location
                samtools sort -@ 3 ${file%.fastq}.unique.bam -o ${file%.fastq}.unique.sorted.bam 2>> ${file%.fastq}.$rand.log
                rm ${file%.fastq}.unique.bam
                
                #more info: https://www.encodeproject.org/data-standards/terms/#library
                #Used for checking library complexity in terms of duplicates dependent on genomic location
                #PBC1 is calculated as number of locations with only one mapping divided by the number of distinct locations with any number of mappings
                #PBC2 is calculated as the number of locations with only one mapping read divided by the number of locations with exactly two reads mapping
                #NRF is calculated as the number of locations with any number of mappings divided by the total number of reads (that weren't multimapped
                mkdir -p qc
                mkdir -p qc/PBC_NRF
                #grab the pbc1, pbc2 and nrf scores before files are rmed forever
                samtools sort -@ 3 -n ${file%.fastq}.unique.sorted.bam | samtools view -bf 0x2 - | bedtools bamtobed -bedpe -i - 2> /dev/null | awk -v OFS="\t" '$1!="chrM" && $4!="chrM"{print $1,$2,$4,$6,$9,$10}' | sort | uniq -c | awk -v OFS="\t" '$1==1{m1++} $1==2{m2++} {md++} {mt+=$1} END {if(mt==0){mdt=-1} else {mdt=md/mt} if(md==0) {m1d=-1} else {m1d=m1/md} if(m2==0) {m12=-1} else {m12=m1/m2} printf "Total Uniquely Mapped Reads:\t%d\nTotal Distinct Genomic Locations:\t%d\nTotal Non-duplicated Reads:\t%d\nTotal 2-plicated reads:\t%d\nNRF:\t%f\nPBC1:\t%f\nPBC2:\t%f",mt,md,m1,m2,mdt,m1d,m12}' > qc/PBC_NRF/${file/.fastq/_PBC_NRF.txt}
                
                #remove dublicates
                samtools rmdup ${file%.fastq}.unique.sorted.bam ${file%.fastq}.bowtie2.bam 2>> ${file%.fastq}.$rand.log
                rm ${file%.fastq}.unique.sorted.bam
                bam=${file%.fastq}.bowtie2.bam
                samtools index $bam


            else
                echo "alignment from" $file "already exists"
        fi
        if [ ! -s $bam ]
            then
                echo "file" $bam "does not exist, whoops"
                exit 1
        fi
        exit 0
fi 

echo "Unknown file extension"
exit 1
