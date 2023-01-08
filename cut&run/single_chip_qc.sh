#! /bin/bash

cd $1
file=$2
genome=$3

if [ ${file: -3} == ".gz" ]
	then
		file=${file%.fastq.gz}.bowtie2.bam
fi

if [ ${file: -6} == ".fastq" ]
	then
		file=${file%.fastq}.bowtie2.bam
fi

if [ ! -s $file ]
	then
		echo "file" $file "does not exist"
		exit 1
fi

if [ ${file: -4} == ".bam" ]
	then
		echo "running qc on" $file
		base=${file%.bam}
		mkdir -p qc
		mkdir -p qc/base-stats
		echo "running bam_stat on" $file
		if [ ! -s qc/base-stats/$base.base-stats.txt ]
			then
			bam_stat.py -i $file > qc/base-stats/$base.base-stats.txt
		fi
		#details: http://rseqc.sourceforge.net/
		#The nucleotide percentage per bp in read
		mkdir -p qc/NVC
        echo "running NVC on" $file
        if [ ! -s qc/NVC/$base.NVC_plot.pdf ]
            then
            read_NVC.py -i $file -o qc/NVC/$base
        fi
        #details http://rseqc.sourceforge.net/
        # Distribution of GC content throughout aligned reads
        mkdir -p qc/GC
        echo "running GC on" $file
        if [ ! -s qc/GC/$base.GC_plot.pdf ]
            then
            read_GC.py -i $file -o qc/GC/$base
        fi
        #details: https://genome.cshlp.org/content/22/9/1813.long
        # Important for ChIP only, not ATAC
        # Essentially, takes the reads on each the Watson and the Crick strands and calculates the Pearson correlation after shifting the Crick strand k base pairs. When plotting correlation vs shift amount, two peaks should appear, one at the read length and one at the fragment length. If the ratio between fragment-length cross-correlation and the background (NSC) is >1.05 and the ratio between the fragment-length peak and read-length peak is >0.8, should be good.
        mkdir -p qc/cross-correlation
        echo "running Cross_Correlation on" $file
        if [ ! -s qc/cross-correlation/$base.log ]
            then
            echo -e "Filename\tnumReads\testFragLen\tcorr_estFragLen\tPhantomPeak\tcorr_phantomPeak\targmin_corr\tmin_corr\tNormalized SCC (NSC)\tRelative SCC (RSC)\tQualityTag" > qc/cross-correlation/$base.log
            Rscript run_spp.R -rf -odir=qc/cross-correlation -p=2 -savp -out=qc/cross-correlation/$base.log -c=$file
        fi
		exit 0
fi 

echo "Unknown file extension"
exit 1
