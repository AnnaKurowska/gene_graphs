##some notes for you Ania
#

#coding_seqs <- readDNAStringSet(orf_fasta_file)--why are we not using fasta file?


setwd("/Users/Ania/Desktop/Szkoła/4th year/Dissertation/gene_graphs/")

source("read_count_functions.R")
source("generate_stats_figs.R") #this doesn't work


library(rhdf5)
library(tidyr)
library(dplyr)
library(magrittr)
library(purrr)
library(ggplot2)

dataset <- "G-Sc_2014"
hd_file <- "Input/WTnone.h5"

# prepare files, opens hdf5 file connection
hdf5file <- rhdf5::H5Fopen(hd_file) # filehandle for the h5 file (first column gives the names and the second gives the type of element listed)

class(hdf5file)
gene_names <- rhdf5::h5ls(hdf5file,
#"H5IdComponent"attr(,"package") "rhdf5"-->how can i deal with that?


# next step is to test with these 5 genes only instead of all genes
#we can use:--> I actually think it's not at this level, because we need to combine both gff and h5 together, that should happen at the later level instead

gene_names <- rhdf5::h5ls(hdf5file, recursive = 1)$name 

# function to get data matrix of read counts for gene and dataset from hdf5file
GetGeneDatamatrix <- function(gene, dataset, hdf5file) {
  hdf5file %>%
    rhdf5::H5Dopen(
      name = paste0("/", gene, "/", dataset, "/reads/data")
    ) %>%
    rhdf5::H5Dread() %>%
    return()
}
#H5Dopen:H5Dopen return an object of class H5IdComponent represinting a H5 dataset identifier. 
#H5DREAD: H5Dread returns an array with the read data.
#MAYBE I SHOULD BE LOOKING AT MODIFYING THIS FUNCTION BC IT OUTPUTS LISTS

#ALSO HERE WE'RE INTRODUCING ALL GENES FROM THE FILE,HOW TO SUSBET IT TO ONLY SELECT GENES("YCR012W","YEL009C","YOR303W","YOL130W","YGR094W") OUT OF 5812 GENES). PERHAPS I COULD DO IT AT AN EARLIER STAGE==GENE NAMES)
allGenesDataMatrix <- lapply(gene_names,
                             function(gene) 
                               GetGeneDatamatrix(gene, dataset, hdf5file)
)



# --- 

orf_gff_file <- "input/yeast_CDS_w_250utrs.gff3"
# WHY ARE WE NOT USING THE FASTA FILE?

# convert GFF to tidy dataframe
readGFFAsDf <- purrr::compose(
  rtracklayer::readGFFAsGRanges,
  data.frame, 
  as_tibble,
  .dir = "forward" # functions called from left to right
)

# read in positions of all exons/genes in GFF format and convert to tibble data frame 
# gives 17436 obs by 10 vars
gff_df <- readGFFAsDf(orf_gff_file)
#SHOULD WE SELECT GFF FOR SPECIFIC GENES AS WELL?

# > gff_df
# # A tibble: 17,436 x 10
# seqnames  start   end width strand source      type  score phase Name     
# <fct>     <int> <int> <int> <fct>  <fct>       <fct> <dbl> <int> <chr>    
# 1 YAL068C       1   250   250 +      rtracklayer UTR5     NA    NA YAL068C  
# 2 YAL068C     251   613   363 +      rtracklayer CDS      NA    NA YAL068C  
# 3 YAL068C     614   863   250 +      rtracklayer UTR3     NA    NA YAL068C


#new<-filter(gff_df, gff_df$seqnames %in% gene_names2)

# extract start locations for each gene from GFF tidy dataframe for CDS only
getCDS5start <- function(name, gffdf, ftype="CDS", fstrand="+") {
  gffdf %>% 
    filter(type==ftype, Name == name, strand == fstrand) %>% 
    pull(start) %>% 
    min 
}

#posn_5start <- GetCDS5start(gene, gff_df) WHY IS THAT COMMENTED?

# function to get matrix of read counts between specific positions
# (from n_buffer before start codon to nnt_gene after start codon)
# for gene and dataset from hd5 file hdf5file, using UTR5 annotations in gff
GetGeneDatamatrix5start <- function(gene, dataset, hdf5file, 
                                    posn_5start, n_buffer, nnt_gene) {
  data_mat_all <- GetGeneDatamatrix(gene, dataset, hdf5file)
  # @ewallace: replace this by gff_df?
  # n_utr5 <- BiocGenerics::width(gff[gff$type == "UTR5" & gff$Name == gene])
  
  # if n_buffer bigger than length n_utr5, pad with zeros:
  if (posn_5start > n_buffer) {
    # if posn_5start bigger than n_buffer
    n_left5 <- posn_5start - n_buffer # column to start from (5'end)
    zeropad5_mat <- matrix(0, nrow = nrow(data_mat_all), ncol = 0)
  } else {
    # if length n_utr5 less than n_buffer
    n_left5 <- 1 # column to start from (5'end)
    zeropad5_mat <- matrix(
      0, 
      nrow = nrow(data_mat_all), 
      ncol = (n_buffer - posn_5start + 1 )
    )
  }
  n_right3 <- posn_5start + nnt_gene - 1 # column to end with (3'end)
  data_mat_5start <- data_mat_all[, n_left5:n_right3]
  return(cbind(zeropad5_mat, data_mat_5start))
}

#sums the list of data matrices
TidyDatamatrix <- function(data_mat, startpos = 1, startlen = 1) {
  # CHECK startpos/off-by-one
  positions <- startpos:(startpos + ncol(data_mat) - 1)
  readlengths <- startlen:(startlen + nrow(data_mat) - 1)
  data_mat %>%
    set_colnames(positions) %>%
    as_tibble() %>%
    mutate(ReadLen = readlengths) %>%
    gather(-ReadLen, key = "Pos", value = "Counts", convert = FALSE) %>%
    mutate(Pos = as.integer(Pos), Counts = as.integer(Counts))
}

# get gene and position specific total counts for all read lengths #i cant identify the genes

gene_poslen_counts_5start_df <-
  lapply(gene_names,
         function(gene) 
           GetGeneDatamatrix5start(gene,
                                   dataset,
                                   hdf5file,
                                   posn_5start = GetCDS5start(gene, gff_df),
                                   n_buffer = 25,
                                   nnt_gene = 50)
  ) %>%
  Reduce("+", .) %>% 
  TidyDatamatrix(startpos = -50 + 1, startlen = 10) 

### attemps to change gene_poslen (matrix) into df
data.frame(t(sapply(gene_poslen_counts_5start_df,c)))
data.frame(Reduce(rbind, gene_poslen_counts_5start_df))
str(enframe(gene_poslen_counts_5start_df))

class(gene_poslen_counts_5start_df) # "tbl_df"     "tbl"        "data.frame"

##

plot_ribogrid <- function(b) {
  ggplot(data = b, aes(x = Pos, y = ReadLen, fill = Counts)) +
    geom_tile() +
    scale_fill_gradient("count", low = "white", high = "darkblue") +
    theme(panel.grid = element_blank()) +
    labs(x = "position of read 5' end", y = "read length")
}


##
#that gives the kind of graph we want to get 
gene_poslen_counts_5start_df %>%
  barplot_ribogrid() %>%
  ggsave(
    filename = file.path(output_dir, paste0(output_prefix, "startcodon_ribogridbar.pdf")),
    width = 6, height = 5
  )


