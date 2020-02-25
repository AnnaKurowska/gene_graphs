setwd("/Users/Ania/Desktop/Szkoła/4th year/Dissertation/gene_graphs/")

library(rhdf5)
library(tidyr)
library(dplyr)
library(magrittr)
library(purrr)

dataset <- "G-Sc_2014"
hd_file <- "G-Sc_2014/output/WTnone/WTnone.h5"

# prepare files, opens hdf5 file connection
hdf5file <- rhdf5::H5Fopen(hd_file) # filehandle for the h5 file


gene_names <- rhdf5::h5ls(hdf5file, recursive = 1)$name
# next step is to test with these 5 genes only instead of all genes
test_orfs<- c("YCR012W","YEL009C","YOR303W","YOL130W","YGR094W")

# function to get data matrix of read counts for gene and dataset from hdf5file
GetGeneDatamatrix <- function(gene, dataset, hdf5file) {
  hdf5file %>%
    rhdf5::H5Dopen(
      name = paste0("/", gene, "/", dataset, "/reads/data")
    ) %>%
    rhdf5::H5Dread() %>%
    return()
}

# --- 

orf_gff_file <- "G-Sc_2014/input/yeast_CDS_w_250utrs.gff3"

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
# > gff_df
# # A tibble: 17,436 x 10
# seqnames  start   end width strand source      type  score phase Name     
# <fct>     <int> <int> <int> <fct>  <fct>       <fct> <dbl> <int> <chr>    
# 1 YAL068C       1   250   250 +      rtracklayer UTR5     NA    NA YAL068C  
# 2 YAL068C     251   613   363 +      rtracklayer CDS      NA    NA YAL068C  
# 3 YAL068C     614   863   250 +      rtracklayer UTR3     NA    NA YAL068C


# extract start locations for each gene from GFF tidy dataframe for CDS only
GetCDS5start <- function(name, gffdf, ftype="CDS", fstrand="+") {
  gffdf %>% 
    dplyr::filter(type==ftype, Name == name, strand == fstrand) %>% 
    dplyr::pull(start) %>%  # pull() pulls out single variable
    min 
}

#posn_5start <- GetCDS5start(gene, gff_df)

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

# get gene and position specific total counts for all read lengths

gene_poslen_counts_5start_df <-
    lapply(test_orfs[1],
           function(gene) 
             GetGeneDatamatrix5start(gene,
                                     dataset,
                                     hdf5file,
                                     posn_5start = GetCDS5start(gene, gff_df),
                                     n_buffer = 25,
                                     nnt_gene = 50)
    ) %>%
  Reduce("+", .) %>% # sums the list of data matrices
  TidyDatamatrix(startpos = -50 + 1, startlen = 10) 

#I've only chosen read lengths between 24 and 30 (y) and then i averaged the read count for each position (x)

y<-gene_poslen_counts_5start_df %>%
  filter(ReadLen >=24 & ReadLen <= 30)

x<-y %>%
  group_by(Pos) %>%
  summarise(avg=mean(Counts))  %>%
  as.data.frame %>%
  arrange(x$Pos)

view(x)
class(x)


    # Runnning these next few lines seems to sum the counts from ALL genes, 
    # and squishes these down into all-gene totals, which isn't what we want.
    # %>%
    # Reduce("+", .) %>% # sums the list of data matrices
    # TidyDatamatrix(startpos = -50 + 1, startlen = 10) 
