#! /bin/bash

cd $1
file=$2
genome=$3
macs14=$4

mkdir -p qc
mkdir -p qc/fastqc
mkdir -p fastq

if [ ${file: -3} == ".gz" ]
    then
        file=${file%.gz}
fi

if [ ${file: -6} == ".fastq" ]
    then
        if [ ! -s $file ]
            then
            echo "file" $file "does not exist"
            exit 1
        fi
        #see https://www.bioinformatics.babraham.ac.uk/projects/fastqc/ for more info
        echo "running fastqc on" $file
        if [ ! -s qc/fastqc/${file%.fastq}_fastqc.zip ]
            then
            fastqc --noextract -j java -o qc/fastqc $file
        fi
        mv $file fastq/
        gzip fastq/$file
        file2=${file/1_val_1/2_val_2}
        if [ ! -s qc/fastqc/${file2%.fastq}_fastqc.zip ]
            then
            fastqc --noextract -j java -o qc/fastqc $file2
        fi
        mv $file2 fastq/
        gzip fastq/$file2
        file=${file%.fastq}.bowtie2.bam
fi

if [ ! -s $file ]
    then
        echo "file" $file "does not exist"
        exit 1
fi

if [ ${file: -4} == ".bam" ]
    then
        base=${file%.bam}
        if [ "$genome" = "mm9" ] || [ "$genome" = "mm10" ]
			then shortgenome="mm"
		fi 
		if [ "$genome" = "hg19" ]
			then shortgenome="hs"
		fi
		mkdir -p qc/macs
		if [ "$macs14" = "yes" ]
			then
				if [ ! -s qc/macs/$base.macs14_peaks.bed ]
					then macs callpeak -t $file -f BAMPE -n qc/macs/$base.macs14 -g $shortgenome
					#Rscript qc/macs/$base.macs14_model.r
				fi
			else
				if [ ! -s qc/macs/$base.macs2_peaks.narrowPeak ]
					then macs2 callpeak -t $file -f BAMPE -n qc/macs/$base.macs2 -q 0.01 -g $shortgenome
					#Rscript qc/macs/$base.macs2_model.r
				fi
		fi
        echo "running bamCoverage on" $file
        fraglen=$(grep "# d =" qc/macs/${base}.macs*_peaks.xls | awk '{print $NF}')
        mkdir -p qc/bw
        if [ ! -s qc/bw/$base.bw ]
            then
            bamCoverage --binSize 25 -p 4 --normalizeUsing RPKM --extendReads $fragLength -b $file -of bigwig -o qc/bw/$base.bw
        fi
        if [ ! -s qc/bw/$base.10x30.bw ]
            then
                case $genome in
                    mm9)
                        bamCoverage --binSize 10 --smoothLength 30 -p 4 --normalizeUsing RPGC --effectiveGenomeSize 2150570000 --extendReads $fragLength -b $file -of bigwig -o qc/bw/$base.10x30.bw
                        ;;
                    mm10)
                        bamCoverage --binSize 10 --smoothLength 30 -p 4 --normalizeUsing RPGC --effectiveGenomeSize 2730871774 --extendReads $fragLength -b $file -of bigwig -o qc/bw/$base.10x30.bw
                        ;;
                    hg19)
                        bamCoverage --binSize 10 --smoothLength 30 -p 4 --normalizeUsing RPGC --effectiveGenomeSize 2451960000 --extendReads $fragLength -b $file -of bigwig -o qc/bw/$base.10x30.bw
                        ;;
                    *)
                        bamCoverage --binSize 10 --smoothLength 30 -p 4 --normalizeUsing RPKM --extendReads $fragLength -b $file -of bigwig -o qc/bw/$base.10x30.bw
                        ;;
                esac
        fi
        exit 0
fi 

echo "Unknown file extension"
exit 1
