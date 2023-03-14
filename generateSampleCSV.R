#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)!=1) {
  stop("Usage: Rscript generateSampleCSV.R <fastq_dir>", call.=FALSE)
}

fastq_dir <- args[1]

# fastq_dir <- "LS180/"

tmp.all <- list.files(fastq_dir)

tmp.fqs <- tmp.all[grep("fq|fastq", tmp.all)]

tmp.format <- ""
if (sum(!grepl("fastq", tmp.fqs)) == 0) {
    tmp.format <- "fastq"
} else if (sum(!grepl("fq", tmp.fqs)) == 0) {
    tmp.format <- "fq"
} else {
    stop("Unknown format!!! file name should be fastq or fq!!!")
}

message(sprintf("Your file format is: %s", tmp.format))

fq1 <- tmp.fqs[grep(paste0("_1.",tmp.format), tmp.fqs)]

fq2 <- gsub(paste0("_1.",tmp.format), paste0("_2.",tmp.format), fq1)

if (sum(!c(fq1, fq2) %in% tmp.fqs) != 0) {
    stop("fq1 and fq2 generation error!!!")
}

sampleName <- unlist(lapply(fq1, function(x) {
    strsplit(x, split = "_")[[1]][1]
}))

SampleCSV <- data.frame(sampleName=sampleName, fq1=paste(getwd(),fastq_dir,fq1,sep = "/"),
                        fq2=paste(getwd(),fastq_dir,fq2,sep = "/"))

write.table(SampleCSV, file = paste0(fastq_dir,".sample.csv"), row.names = F, col.names = F, quote = F, sep = ",")
