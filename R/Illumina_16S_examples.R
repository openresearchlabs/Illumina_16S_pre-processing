

## Illumina 16S pre-processing
## Usage examples 
## 20 April 2015


## ===================================================================================
## Preparations
## Assumes raw data does not contain adapter or barcode sequences
## FWD reads are indicated fy "_R1_" in the filename, reverse by "_R2_"
## ===================================================================================

# set working directory to where reads data is
setwd('$PWD/Illumina_preprocessing/data')

# source functions
source('$PWD/Illumina_preprocessing/Illumina_16S_functions.R')

# generate sample names
sample.names <- sapply(list.files(pattern='_R1_'), function(x) strsplit(x, '-', fixed=T)[[1]][1])


## ===================================================================================
## Merge paired reads
## FWD reads are indicated fy "_R1_" in the filename, reverse by "_R2_"
## Writing nonmerged FWD and REV reads to fasta and fastq is set by output.notmerged=T
## ===================================================================================

mergeReads(forward = list.files(pattern='_R1_.*fastq$'), # list forward fastq files
           reverse = list.files(pattern='_R2_.*fastq$'), # list reverse fastq files
           out.names = sample.names,
           min.merged.length = 450,
           min.overlap = 30,
           max.mismatches = 10,
           output.notmerged = FALSE) # set this to TRUE to output nonmerged reads in fasta and fastq format

# quality stats for merged and fwd/rev reads
runEEStats('sample001_merged.fastq')
runEEStats('sample001-ATATGGGA-exampledata_S117_L001_R1_001.fastq')
runEEStats('sample001-ATATGGGA-exampledata_S117_L001_R2_001.fastq')


## ===================================================================================
## Quality filtering 
## ===================================================================================

filterReads(files = list.files(pattern='_merged.fastq$'), # list merged fastq files
            trucation.quality = 10,
            max.expected.error = 1.0,
            min.length = 450, 
            out.names = sample.names)

# running stats for a filtered file
runEEStats('sample002_filtered.fastq')

# plot histogram of read lengths for a filtered file
hist(width(readDNAStringSet('sample001_filtered.fasta')))



## ===================================================================================
## Quality filtering for non-merged reads
## 
## ===================================================================================

filterReads(files = list.files(pattern='_001.fastq$'), # list fwd and rev fastq files
            trucation.quality = 10,
            max.expected.error = 1.0,
            min.length = 200, 
            out.names = gsub('_001', '_001_filtered', list.files(pattern='_001.fastq$'))) # output filenames generated from original raw fastq filenames

# Combining fwd and rev to the same file, make reverse complement of rev
combineSE(files = list.files(pattern='*.R.*_filtered.fasta$'), # lists filtered fwd and rev fastq files
          out.names=sample.names)


## ===================================================================================
## Chimera filtering
## ===================================================================================

chimeraFilter(files = list.files(pattern='_filtered.fasta$'), # lists filtered FASTA files
              out.names = sample.names,
              ref.db = '/Volumes/MyBook/jarmo_projects/DBs/rdp_gold.fasta')


## ===================================================================================
## Equalise read numbers (not necessarily needed!)
## ===================================================================================

plotReadNumbers(files = list.files(pattern='_nonchimeric.fasta$'), # lists chimera filtered fasta files
                out.names = sample.names)
  

equalizeReadNumbers(files = list.files(pattern='_nonchimeric.fasta$'), # same as above
                    cutoff = 3614, # equalise to this number
                    out.names = sample.names) 


## ===================================================================================
## Generate files for Qiime
## Input equalised or non.equalised filtered fasta files
## ===================================================================================

combineSamples(files = list.files(pattern='_nonchimeric.fasta$'), # lists nonchimeric fasta files
               out.names = sample.names,
               out.fasta.name = 'all_noneq') # output filename for non-equalised fasta





