#!/bin/bash
# POSIX

directions (){
	echo ""
	echo "Chip analysis pipeline by Nick - Similar to mapmaker"
	echo "Counts BAM reads over set regions and creates associated plots"
	echo "Usage:"
	echo "makemapper -g mm9 -b ROI.bed -k -q -c 10 -r -p 90 -m 2 -n 2 -o makemapper_out -s promoters.bed -i file1.bam file2.bam file3.bam"
	echo "-g | --genome:      REQUIRED: genome"
	echo "-b | --bed:         OPTIONAL: ROI file on which to count, macs2 peaks are merged by default"
	echo "-s | --subtract:    OPTIONAL: Region to subract from ROI (i.e. promoters)"
	echo "-k | --rpkm:        Run RPKM normalization instead of RPM"
	echo "-q | --quantile:    Run quantile normalization in R"
	echo "-c | --clusters:    Number of kmeans clusters, 1 by default"
	echo "-r | --nohier:      Disable hierarchical clustering sort of samples"
	echo "-p | --percentcut:  Percentage of least differential regions to cut. Default 90%"
	echo "-m | --minval:      Minimum RPM/RPKM in sample to keep region, combine with -n"
	echo "-n | --minnum:      Minimum number of samples reaching minval to keep region, combine with -m"
	echo "-o | --out:         Output name for pairwise comparisons and pdf"
	echo "-i | --bams:        REQUIRED: Bam files to count over regions"
	echo "## Will not overwrite created files, quick to run after a single change"
}

#set arguments to defaults or blanks
genome=""
bed=""
rpkm="false"
quant=0
perc=90
minval=-1
minnum=2
outname="makemapper_out"
subbed=""
hier=1
kclusts=1
bams=""



#parse arguments
while :; do
	case $1 in
		-h|--help)
			directions
			exit
			;;
                -d|--directory)
			cd $2
                        shift
			;;
		-g|--genome)
			genome=$2
			shift
			;;
		-i|--bams)
			#first bam file
			if  [[ ${2} != "-"* ]] ;
			then
				bams=$2
			else
				printf '-i requires bam files.\n' >&2
				exit 1
			fi
			shift
			#the other bam files
			while [[ ${2} != "-"* ]] && [ ! -z "$2" ]
			do
				bams=$bams" "$2
				shift
			done
			;;
		-b|--bed)
			if [ -n "$2" ]
				then bed=$2
				shift
			else
				printf '-b requires an existing bed file.\n' >&2
				exit 1
			fi
			;;
		-k|--rpkm)
			rpkm="true"
			;;
		-q|--quantile)
			quant=1
			;;
		-c|--clusters)
			kclusts=$2
			shift
			;;
		-r|--nohier)
			hier=0
			;;
		-p|--percentcut)
			perc=$2
			shift
			;;
		-m|--minval)
			minval=$2
			shift
			;;
		-n|--minnum)
			minnum=$2
			shift
			;;
		-o|--out)
			outname=$2
			shift
			;;
		-s|--subtract)
			if [ -n "$2" ]
				then subbed=$2
				shift
			else
				printf '-s requires an existing bed file.\n' >&2
				exit 1
			fi
			;;
		--)
			shift
			break
			;;
		-?*)
			printf 'unknown option '$1'\n' >&2
			directions
			exit
			;;
		*)
			break
			;;
	esac
	shift
done

#check if genome is set to acceptable value
if [ "$genome" != "mm9" ] && [ "$genome" != "hg19" ] && [ "$genome" != "mm10" ]
	then
	echo ""
	echo "Missing required argument -g [mm9,mm10,hg19]"
	echo ""
	directions
	exit
fi

#main loop to create bed regions
#if the bed argument is empty....
if [ ! -n "$bed" ]
	then 
	#set bed argument to use this file from now on
	bed=${outname/.bam/}".merge.bed"
	#loop to run macs on all bam files, with appropriate genome words
	for var in $bams
		do temp=${var%.bam}
		if [ "$genome" = "mm9" ] || [ "$genome" = "mm10" ]
			then shortgenome="mm"
		fi 
		if [ "$genome" = "hg19" ]
			then shortgenome="hs"
		fi 
		echo $temp
		if [ ! -s "qc/macs/$(basename $temp).macs2_peaks.narrowPeak" ]
			then /apps/python-2.7.9/bin/macs2 callpeak -t $var -f BAM -n qc/macs/$(basename $temp).macs2 -q 0.01 -g $shortgenome
			/ifs/labs/shivdasani/sv681/ref/R-3.4.0/bin/Rscript qc/macs/$(basename $temp).macs2_model.r
		fi
	done
	#create the bed file of interest by merging MACS results and sorting.  Subtract a bed file (promoter?) if the option is set
	#The special sort file (set by faidx) is NOT normal bed sorting, it is the order used in most BAM files and is necessary for bedtools intersect -c later
	if [ ! -n "$subbed" ]
		then cat qc/macs/*.narrowPeak | awk '{print $1,$2,$3}' OFS="\t" | /apps/bedtools2-2.25.0/bin/bedtools sort -faidx /ifs/labs/shivdasani/sv681/ref/chr_bad_order.txt -i - | /apps/bedtools2-2.25.0/bin/bedtools merge -i - > $outname".merge.bed"
	else
		cat qc/macs/*.narrowPeak | awk '{print $1,$2,$3}' OFS="\t" | /apps/bedtools2-2.25.0/bin/bedtools sort -faidx /ifs/labs/shivdasani/sv681/ref/chr_bad_order.txt -i - | /apps/bedtools2-2.25.0/bin/bedtools subtract -A -a - -b $subbed | /apps/bedtools2-2.25.0/bin/bedtools merge -i - > $outname".merge.bed"
	fi
else
#if the bed file was given, sort (according to BAM order) and subtract bed file if the option is set
	if [ ! -n "$subbed" ]
	then
		/apps/bedtools2-2.25.0/bin/bedtools sort -faidx /ifs/labs/shivdasani/sv681/ref/chr_bad_order.txt -i $bed > $bed.sort.bed
		bed=$bed".sort.bed"
	else
		/apps/bedtools2-2.25.0/bin/bedtools sort -faidx /ifs/labs/shivdasani/sv681/ref/chr_bad_order.txt -i $bed | /apps/bedtools2-2.25.0/bin/bedtools subtract -A -a - -b $subbed > $bed.sort-limit.bed
		bed=$bed".sort-limit.bed"
	fi
fi



#index all of the bam files
for var in $bams
	do if [ ! -s $var.bai ]
		then samtools index $var
	fi
done

#create the original count matrix on the bed file of interest with BAM reads
if [ ! -s $outname.rpm-counts-matrix.tsv ]
	#we will be modifying temp with additional columns for each bam file
	then cp $bed temp
	#prime the header by predicting column 1 as 'locus'
	printf "locus" > header.txt
	for var in $bams
		#continue the header construction per column as columns are added to the BED file
		do printf "\t"$(basename $var) >> header.txt
		#grab the size of the current BAM file, for RPM normalization in a moment
		normval=$(samtools view -c $var)
		#add a column to the bed file (temp) with BAM read information, and then divide this 'last column' by the number of reads to get RPM
		/apps/bedtools2-2.25.0/bin/bedtools intersect -c -g /ifs/labs/shivdasani/sv681/ref/chr_bad_order.txt -sorted -a temp -b $var | awk -v normval=$normval '{$NF=$NF*1000000/normval;print}' OFS="\t" > temp.2
		mv temp.2 temp
	done
	#finalize the header
	printf "\n" >> header.txt
	#condence the locus into a single column (remove first two tabs), add the header, remove chrY and chrM and chrRANDOMS, matrix done!
	sed 's/\t/-/' temp | sed 's/\t/-/' | cat header.txt - | grep -v 'chrM' | grep -v 'chrY' | grep -v 'random' > $outname.rpm-counts-matrix.tsv
	rm header.txt
	rm temp
fi


#if RPKM is desired, create a new matrix by dividing EACH column by the length of the locus and multiplying by 1000.  Have to unpack and repack the 'matrix' format
if [ ! -s $outname.rpkm-counts-matrix.tsv ] && $rpkm
	then awk 'NR==1{print}' $outname.rpm-counts-matrix.tsv > header.txt
	awk 'NR!=1{val=1000/($3-$2); for (i=4;i<NF+1;i++) $i=$i*val; print}' FS="[-\t ]" OFS="\t" $outname.rpm-counts-matrix.tsv | sed 's/\t/-/' | sed 's/\t/-/' | cat header.txt - > $outname.rpkm-counts-matrix.tsv
	rm header.txt
fi

#limit the regions to those that have a certain value (minval) at least n (minnum) number of times, and run the R script with all relevant arguments. Input file depends on RPKM argument
if $rpkm
	then awk -v minnum=$minnum -v minval=$minval '{num=0;for (i=2;i<NF+1;i=i+1) if ($i > minval) num=num+1} num>minnum-1{print}' $outname.rpkm-counts-matrix.tsv > $outname.rpkm-limit.tsv
	/ifs/labs/shivdasani/sv681/ref/R-3.4.0/bin/Rscript /ifs/labs/shivdasani/sv681/ref/bin/mapmaker.r $outname $outname.rpkm-limit.tsv $quant $perc $hier $kclusts
else
	awk -v minnum=$minnum -v minval=$minval '{num=0;for (i=2;i<NF+1;i=i+1) if ($i > minval) num=num+1} num>minnum-1{print}' $outname.rpm-counts-matrix.tsv > $outname.rpm-limit.tsv
	/ifs/labs/shivdasani/sv681/ref/R-3.4.0/bin/Rscript /ifs/labs/shivdasani/sv681/ref/bin/mapmaker.r $outname $outname.rpm-limit.tsv $quant $perc $hier $kclusts
fi

#this file just added more mess,though you may want to keep it
rm *limit.tsv
















