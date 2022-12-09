#! /bin/bash

directions (){
	echo ""
	echo "ChIP pipeline by Nick, modified by Shariq"
	echo "Runs bowtie2 and supplemental programs"
	echo "Usage:"
	echo "chip_core_PE.sh -q -t -g [mm9,mm10,hg19] sample1 sample2 sample3"
    echo "-h: display this helpful information"
	echo "-g: REQUIRED - genome"
	echo "-q: Run QC on bam"
    echo "-m: Run makemapper on all bams"
    echo "-t: Run trimming on fastqs"
	echo "## Will NOT overwrite bam files that already exist and are not empty"
}

module load bowtie2-2.3.5-gcc-5.4.0-cxapxq6 samtools-1.9-gcc-5.4.0-jjq5nua fastqc-0.11.7-gcc-5.4.0-hjgfsw2 bedtools2-2.27.1-gcc-5.4.0-iouh4nk

qc=0
misc=0
genome=0
mm=0
macs14="no"
trim="no"

while getopts "g:hmq1t" OPTION
do
	case $OPTION in
		h)
			directions
			exit 0
			;;
		q)
			qc=1
			;;
        m)
            mm=1
            ;;
		g)
			genome=$OPTARG
			;;
        1)
			macs14="yes"
			;;
        t)
            trim="yes"
            ;;
	esac
done



#check required args
if [ "$genome" != "bt6" ] && [ "$genome" != "mm9" ] && [ "$genome" != "hg19" ] && [ "$genome" != "mm10" ]
	then
	echo ""
	echo "Missing required argument -g [mm9,mm10,hg19]"
	directions
	exit
fi

#create random id for the run
rand=$RANDOM

#id to wait for job on all output
waitall=""
waitalign=""

#initialize list of bams
bams=""

#parse remaining args (the files and groups)
theArgs=$(echo $@ | awk -v optind=$OPTIND '{for (i=optind; i<NF+1; i++) printf $i" "}' | awk '$1!="-s"{print}' RS="[ \n]")

#check that files exist
for line in $theArgs; do
	if [ ! -s $line ]
		then
		    echo ""
		    echo "File: $line not found"
		    directions
		    exit
	fi
done

#make log directory
mkdir -p logs

#Loop to analyze all files seperately
for line in $theArgs; do
	if [ ${line: -3} == ".gz" ] || [ ${line: -6} == ".fastq" ]
		then
	    	waitone=($(qsub -cwd -v PATH -o ${PWD}/logs -e ${PWD}/logs single_chip_alignment_PE.sh $PWD $line $genome $rand $trim | awk '{print $3","}'));
		    if [ $trim == "yes" ]
                then
                    line=${line/.fastq/_val_1.fastq}
            fi
            bams=$bams" "${line%.gz}
            bams=${bams%.fastq}.bowtie2.bam
		    waitall="$waitall""$waitone";
            waitalign="$waitalign""$waitone";
		    if [ $qc == "1" ]
			    then
			    echo ""
			        waitwon=($(qsub -cwd -v PATH -o ${PWD}/logs -e ${PWD}/logs -hold_jid $waitone single_chip_qc.sh $PWD $line $genome | awk '{print $3","}'))
			        waitall="$waitall""$waitwon";
		    fi
		    if [ $qc == "1" ]
			    then
			        echo ""
			        qsub -cwd -v PATH -o ${PWD}/logs -e ${PWD}/logs -hold_jid $waitone single_chip_misc_PE.sh $PWD $line $genome
		    fi
	fi
done

#run makemapper on all alignments
if [ $mm == "1" ]
        then
            qsub -cwd -v PATH -o ${PWD}/logs -e ${PWD}/logs -hold_jid $waitalign /ifs/labs/shivdasani/sv681/ref/bin/makemapper.sh --directory $PWD -g $genome -s "/ifs/labs/shivdasani/sv681/ref/"$genome"/"$genome".promoter-1k2k.bed" -q -c 10 -i $bams
fi

#run group qc
if [ $qc == "1" ]
	then
	echo "";
	#echo $PWD
	#echo $genome
	#echo $rand
	#echo $theArgs	
	qsub -cwd -v PATH -o ${PWD}/logs -e ${PWD}/logs -hold_jid $waitall chip_group_qc.sh $PWD $genome $rand $theArgs
fi
