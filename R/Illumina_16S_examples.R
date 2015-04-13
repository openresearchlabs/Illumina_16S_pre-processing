

## Illumina 16S pre-processing
## Usage examples 
## 13 April 2015


## ================================================================
## Preparations
## Assumes raw data does not contain adapter or barcode sequences
## ================================================================

# set working dir
setwd('~/Illumina_preprocessing/data')

# source functions
source('~/Illumina_preprocessing/Illumina_16S_functions.R')

# generate sample names
sample.names <- sapply(list.files(pattern='_R1_'), function(x) strsplit(x, '-', fixed=T)[[1]][1])


## ================================================================
## Merge paired reads
## ================================================================

mergeReads(forward = list.files(pattern='_R1_.*fastq$'),
           reverse = list.files(pattern='_R2_.*fastq$'),
           out.names = sample.names,
           min.merged.length = 450,
           min.overlap = 30,
           max.mismatches = 10)


runEEStats('sample001_merged.fastq')
runEEStats('sample001-ATATGGGA-exampledata_S117_L001_R1_001.fastq')
runEEStats('sample001-ATATGGGA-exampledata_S117_L001_R2_001.fastq')


## ================================================================
## Quality filtering of merged reads
## ================================================================

filterReads(files = list.files(pattern='_merged.fastq$'),
            trucation.quality = 10,
            max.expected.error = 1.0,
            min.length = 450, 
            out.names = sample.names)

runEEStats('sample002_filtered.fastq')

hist(width(readDNAStringSet('sample001_filtered.fasta')))



## ================================================================
## Quality filtering for non-merged files
## Combining pairs to the same file
## ================================================================

filterReads(files = list.files(pattern='_001.fastq$'),
            trucation.quality = 10,
            max.expected.error = 1.0,
            min.length = 200, 
            out.names = gsub('_001', '_001_filtered', list.files(pattern='_001.fastq$')))

combineSE(files = list.files(pattern='*.R.*_filtered.fasta$'), out.names=sample.names)


## ================================================================
## Chimera filtering
## ================================================================

chimeraFilter(files = list.files(pattern='_filtered.fasta$'),
              out.names = sample.names,
              ref.db = '/Volumes/MyBook/jarmo_projects/DBs/rdp_gold.fasta')


## ================================================================
## Equalise read numbers (not necessarily needed!)
## ================================================================

plotReadNumbers(files = list.files(pattern='_nonchimeric.fasta$'),
                out.names = sample.names)
  

equalizeReadNumbers(files = list.files(pattern='_nonchimeric.fasta$'),
                    cutoff = 3614,
                    out.names = sample.names) 


## ================================================================
## Generate files for Qiime
## ================================================================

combineSamples(files = list.files(pattern='_nonchimeric.fasta$'),
               out.names = sample.names,
               out.fasta.name = 'all_noneq')





