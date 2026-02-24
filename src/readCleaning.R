# This script handles
# 1. Deduplication
# 2. Filtering of reads containing CA
# 3. UMI & CA removal

library(Biostrings)
library(data.table)

trimmedDir <- snakemake@params[["trimmed_dir"]]
outDir <- snakemake@params[["out_dir"]]

if(!dir.exists(outDir)) dir.create(outDir)

file <- snakemake@input[["trimmed_fastq"]]
unzip(file)
readSeq <- readDNAStringSet(file, format="fastq")
umiDt <- data.table(seq=as.character(readSeq),
                    umi=substr(as.character(readSeq), 1, 8),
                    ca=substr(as.character(readSeq), 9,10))
umiDt[,ind:=1:nrow(umiDt)]
nTot <- nrow(umiDt)
statsFileDt <- data.table(nTot=nTot)
umiDt <- subset(umiDt, ca=="CA")
statsFileDt[,hasCA:=nrow(umiDt)/nTot]
  
umiDt <- unique(umiDt, by="seq")
statsFileDt[,uniqueUMI:=nrow(umiDt)/nTot]
statsFileDt$file <- file

readClean <- readSeq[umiDt$ind]
readClean <- subseq(readClean, start=11)

outFile <- snakemake@output[["filtered_fastq"]]
writeXStringSet(readClean, outFile, format="fastq")
write.table(statsFileDt, file=file.path(dirname(outFile), 
                                        paste("filtering_stats", 
                                              gsub("fastq", "tsv", 
                                              basename(outFile)), 
                                              sep="_")),
            row.names=FALSE, quote=FALSE, sep="\t")