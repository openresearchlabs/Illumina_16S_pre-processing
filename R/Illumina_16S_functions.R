

## ================================================================
## Illumina pre-processing functions
## April 20 2015
## ================================================================


library(Biostrings)


## hardcoded path to Usearch
my.usearch.path <- '$PWD/usearch8.0.1517'


##
## Merge read pairs.
##
mergeReads <- function(forward = list.files(pattern='_R1_.*fastq$'), # list all fastq files containing 'R1'
                       reverse = list.files(pattern='_R2_.*fastq$'), #  
                       out.names, # sample names
                       min.merged.length = 450,
                       min.overlap = 10,
                       max.mismatches = 3,
                       output.notmerged = F, # set this to TRUE to output nonmerged reads in fasta and fastq formats
                       usearch.path = my.usearch.path) {
  
  # call usearch to merge fastq pairs
  if(output.notmerged) {
    for(i in seq_along(forward)) {
      system(paste(usearch.path, ' -fastq_mergepairs ', forward[i], ' -reverse ', reverse[i], ' -fastq_minmergelen ', min.merged.length, 
                   ' -fastq_minovlen ', min.overlap, ' -fastq_maxdiffs ', max.mismatches, ' -fastqout ', out.names[i], '_merged.fastq', 
                   ' -fastqout_notmerged_fwd ', out.names[i], '_notmerged_fwd.fastq', ' -fastqout_notmerged_rev ', out.names[i], '_notmerged_rev.fastq', 
                   ' -fastaout_notmerged_fwd ', out.names[i], '_notmerged_fwd.fasta', ' -fastaout_notmerged_rev ', out.names[i], '_notmerged_rev.fasta', sep=''))
    }
  } else for(i in seq_along(forward)) {
    system(paste(usearch.path, ' -fastq_mergepairs ', forward[i], ' -reverse ', reverse[i], ' -fastq_minmergelen ', min.merged.length, 
                 ' -fastq_minovlen ', min.overlap, ' -fastq_maxdiffs ', max.mismatches, ' -fastqout ', out.names[i], '_merged.fastq', sep=''))
  }
}


##
## Filtering reads based on quality metrics.
##
filterReads <- function(files = list.files(pattern='_merged.fastq$'), # list merged files
                        trucation.quality = 8,
                        max.expected.error = 0.5,
                        min.length = 450,
                        #truncation.length = 1,  # delete after N'th base
                        #stripleft = 1, # delete first N bases
                        out.names, # sample names
                        usearch.path = my.usearch.path ) {
  
  # run fastq quality filtering for paired
  for(i in seq_along(files)) {
    system(paste(usearch.path, ' -fastq_filter ', files[i], ' -fastq_truncqual ', trucation.quality, 
                 ' -fastq_maxee ', max.expected.error, ' -fastq_minlen ', min.length,
                 ' -fastqout ', out.names[i], '_filtered.fastq', 
                # ' -fastq_trunclen ', truncation.length, ' -fastq_stripleft ', stripleft,
                 ' -fastaout ', out.names[i], '_filtered.fasta', sep=''))
  }
} 


##
## Combine *non-merged* paired fasta reads to the same file.
##
combineSE <- function(files = list.files(pattern='_filtered.fasta$'), out.names) {
  fwd <- files[which(grepl('R1', files))]
  rvs <- files[which(grepl('R2', files))]
  if(length(rvs) == 0) cat('\n No input files found \n')
  for(i in 1:length(rvs)) {
    writeXStringSet(c(readDNAStringSet(fwd[i]), reverseComplement(readDNAStringSet(rvs[i]))), 
                    paste(out.names[i], '_R1R2.fasta', sep=''))
  }
}


##
## Check quality stats.
##
runEEStats <- function(file,
                       usearch.path = my.usearch.path) {
  
  # call usearch ee stats function
  sample.name <- gsub('.fastq', '', file)
  system(paste(usearch.path, ' -fastq_eestats ', file, ' -output ', sample.name, '_eestats.txt', sep=''))
  
  # plot EE stats results
  op <- par(mfrow = c(1,2))
  tmp.stats <- read.delim(paste(sample.name, '_eestats.txt', sep=''), header=T, stringsAsFactors=F)
  plot((tmp.stats[, 'Mean_EE']), type='l', lwd=2, col='black', lty='dashed', main='Cumulative Expected Error by position', ylab='Expected Error', xlab='Position')
  lines(tmp.stats[, 'Med_EE'], type='l', lwd=2, col='grey', lty='dashed')
  lines(tmp.stats[, 'Low_EE'], type='l', lwd=2, col='orange', lty='dashed')
  lines(tmp.stats[, 'Min_EE'], type='l', lwd=2, col='red', lty='dashed')
  lines(tmp.stats[, 'Hi_EE'], type='l',  lwd=2, col='forestgreen', lty='dashed')
  legend('topleft', lwd=2, lty='dashed', col=c('black', 'grey', 'orange', 'red', 'forestgreen'), legend=c('Mean EE', 'Median EE', 'Lower quartile EE', 'Min EE', 'Upper quartile EE'))
  plot(tmp.stats[, 'Mean_Q'], col='blue', lwd=2, type='l', main='Mean Phred Q-score by position', ylab='Phred Q-score', xlab='Position')
  par(op)
  
}


##
## Reference based filtering of probable chimeric reads.
##
chimeraFilter <- function(files = list.files(pattern='_filtered.fasta$'),
                          ref.db = '$PWD/rdp_gold.fasta',
                          out.names,
                          usearch.path = my.usearch.path) {
  
  # call usearch to run reference based filtering
  for(i in seq_along(files)) {
    system(paste(usearch.path, ' -uchime_ref ', files[i], ' -db ' , ref.db, ' -strand plus -nonchimeras ',
               paste(out.names[i], '_nonchimeric.fasta', sep=''))) 
  }

}
  

##
## Plotting and printing read numbers per sample.
##
plotReadNumbers <- function(files = list.files(pattern='_nonchimeric.fasta$'),
                            out.names) {
  
  # barplot and table of read numbers in each sample
  nreads <- sapply(files, function(x) length(readDNAStringSet(x))) # count read numbers
  op <- par(mfcol=c(1,1), mfrow=c(1,1))
  barplot(nreads, names.arg=out.names, las=2) # plot read numbers
  par(op)
  data.frame(cbind(out.names, nreads))

}
  

##
## Equalises read numbers to the value specified by 'cutoff'.
## If read number < cutoff, all reads are included.
##
equalizeReadNumbers <- function(files = list.files(pattern='_nonchimeric.fasta$'),
                                cutoff,
                                out.names) {  
  
  # sample randomly the number of reads indicated by 'cutoff'
  for(i in seq_along(files)) {
    tt <- readDNAStringSet(files[i])
    if(length(tt) >= cutoff) {
      keep <- sample(1:length(tt), cutoff)
      writeXStringSet(tt[keep], paste(out.names[i], '_eq.fasta', sep=''))
    } else {
      writeXStringSet(tt, paste(out.names[i], '_eq.fasta', sep=''))
    }  
  }

}
  

## 
## Combining all samples into the same file
## New barcodes are generated for each sample to prevent barcode overlapping issues.
## Generates a Qiime formatted (minimal) mapping file. LinkerPrimer is not intended for use.
## Warns if output files of the same name already exists because data is appended to the files. 
##
combineSamples <- function(files = list.files(pattern='_nonchimeric.fasta$'),
                           out.names,
                           out.fasta.name) {
  
  # output fasta filename
  outfsa <- paste(out.fasta.name, '.fasta', sep='')
  
  # if a combined file exists, notify
  if(length(list.files(pattern=outfsa)) > 0) cat('\n  Warning: a combined fasta file already exists!  \n' )
  
  # generate pseudobarcodes
  new.bc <- unique(sapply(1:2000, function(x) paste(sample(c('A', 'C', 'T', 'G'), 8, replace=T), collapse='')))[1:length(files)]
  names(new.bc) <- out.names

  # merge pseudobarcodes to sequences and combine all samples
  for(i in seq_along(files)) {
    tmp <- readDNAStringSet(files[i])
    tmp <- xscat(DNAStringSet(new.bc[i]), tmp) # merge barcodes
    writeXStringSet(tmp, outfsa, append=T)
  }
  
  # rename sequences
  allpros <- readDNAStringSet(outfsa)
  names(allpros) <- 1:length(allpros)
  writeXStringSet(allpros, outfsa)
  rm(allpros); gc()
  
  # generate Qiime mapping file
  if(length(list.files(pattern='Qiime_mapping.txt$')) > 0) cat('\n  Warning: a Qiime mapping file already exists and will be deleted!  \n' )
  if(length(list.files(pattern='Qiime_mapping.txt$')) > 0) system('rm Qiime_mapping.txt')  
  cat('#SampleID', '\t',  'BarcodeSequence', '\t',  'LinkerPrimerSequence', '\t',  'Treatment', '\t',  'DOB', '\t',  'Description', '\n', file='Qiime_mapping.txt', sep='')
  for(i in seq_along(sample.names)) {
    cat(out.names[i], '\t', new.bc[i], '\t', 'YATGCTGCCTCCCGTAGGAGT', '\t', 'Treatment', '\t', format(Sys.time(), "%d%m%Y"), '\t', 'Description', '\n', append=T, sep='', file='Qiime_mapping.txt')
  }
  
}




