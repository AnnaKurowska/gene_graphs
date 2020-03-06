<<<<<<< HEAD
rm(list=ls())
if (!require(xlsx)){
  install.packages("xlsx")
  library(xlsx)
}

#Planning:

#basically the GetGeneDataMatrix seems fine:
=======
###Set up:

# need this
library(tibble)

nnt_gene<- 50
setwd("/Users/Ania/Desktop/Szkoła/4th year/Dissertation/gene_graphs/")

dataset <- "G-Sc_2014"
hd_file <- "G-Sc_2014/output/WTnone/WTnone.h5"

# prepare files, opens hdf5 file connection
hdf5file <- rhdf5::H5Fopen(hd_file) # filehandle for the h5 file
# formal class H5IDComponent

gene_names <- rhdf5::h5ls(hdf5file, recursive = 1)$name
#  chr [1:5812] "Q0045" "Q0050"

# next step is to test with these 5 genes only instead of all genes
test_orfs<- c("YCR012W","YEL009C","YOR303W","YOL130W","YGR094W")

orf_gff_file <- "G-Sc_2014/input/yeast_CDS_w_250utrs.gff3"
# chr [1:5] "YCR012W" "YEL009C"

# FLIC: moved function up 
readGFFAsDf <- purrr::compose(
  rtracklayer::readGFFAsGRanges,
  data.frame, 
  as_tibble,
  .dir = "forward" # functions called from left to right
)

gff_df <- readGFFAsDf(orf_gff_file)
# > str(gff_df)
# Classes ‘tbl_df’, ‘tbl’ and 'data.frame':	17436 obs. of  10 variables:
#   $ seqnames: Factor w/ 5812 levels "YAL001C","YAL002W",..: 68 68 68 67 67 67 66 66 66 65 ...
# $ start   : int  1 251 614 1 251 479 1 251 2033 1 ...
# $ end     : int  250 613 863 250 478 728 250 2032 2282 250 ...
# $ width   : int  250 363 250 250 228 250 250 1782 250 250 ...
# $ strand  : Factor w/ 3 levels "+","-","*": 1 1 1 1 1 1 1 1 1 1 ...
# $ source  : Factor w/ 1 level "rtracklayer": 1 1 1 1 1 1 1 1 1 1 ...
# $ type    : Factor w/ 3 levels "UTR5","CDS","UTR3": 1 2 3 1 2 3 1 2 3 1 ...
# $ score   : num  NA NA NA NA NA NA NA NA NA NA ...
# $ phase   : int  NA NA NA NA NA NA NA NA NA NA ...
# $ Name    : chr  "YAL068C" "YAL068C" "YAL068C" "YAL067W-A" ...

# used: purrr::compose(rtracklayer::readGFFAsGRanges, data.frame, as_tibble, .dir = "forward")
# FLIC: use this to as_tibble stuff later?

>>>>>>> 1aaa8cea6ea77b8a48e506e09048f49e5ac3be60
GetGeneDatamatrix <- function(gene, dataset, hdf5file) {
  hdf5file %>%
    rhdf5::H5Dopen(
      name = paste0("/", gene, "/", dataset, "/reads/data")
    ) %>%
    rhdf5::H5Dread() %>%
    return()
}

# gff_df %>%
#   filter(type=="UTR5", strand == "+")
# # A tibble: 5,812 x 10
# seqnames  start   end width strand source      type  score phase Name     
# <fct>     <int> <int> <int> <fct>  <fct>       <fct> <dbl> <int> <chr>    
#   1 YAL068C     1   250   250 +      rtracklayer UTR5     NA    NA YAL068C  
# 2 YAL067W-A     1   250   250 +      rtracklayer UTR5     NA    NA YAL067W-A

# mutate() adds new variables, leaves old ones
# transmute() adds new variables, drops old ones

# THIS NEEDS RENAMING :D
testfunction <- function(x){
  x %>%
    filter(type=="UTR5", strand == "+") %>%
    # # A tibble: 5,812 x 10
    mutate(utr5start = (end * -1), utr5end = (start -2))  # FLIC: should this be -2 or -1?
}

# apply the NEWLY RENAMED function ^^^
utr5_data <- testfunction(gff_df)
# > utr5_data
# # A tibble: 5,812 x 4
# Name      utr5start utr5end utr5len
# <chr>         <dbl>   <dbl>   <dbl>
# 1 YAL068C        -250      -1     250
# 2 YAL067W-A      -250      -1     250
# 3 YAL067C        -250      -1     250


# data_mat_all <- GetGeneDatamatrix(test_orfs[1], dataset, hdf5file)
# # > str(data_mat_all)
# # int [1:41, 1:1751] 0 0 0 0

# needs comments
UTR5_length <- function(gene, x) {
  x %>%
    dplyr::filter(Name == gene) %>%
    dplyr::pull(utr5len)
}

# needs comments
Get5UTRstart <- function(gene, x) {
  x %>% 
    dplyr::filter(Name == gene) %>% 
    dplyr::pull(start)
}

# needs comments
UTR5_length <- function(gene, x) {
  x %>% 
    dplyr::filter(Name == gene) %>% 
    dplyr::pull(width)
} 

# needs comments
UTR5_5start <- function(gene, x) {
  x %>% 
    dplyr::filter(Name == gene) %>% 
    dplyr::pull(utr5start)
} 

# needs comments
UTR5_5end <- function(gene, x) {
  x %>% 
    dplyr::filter(Name == gene) %>% 
    dplyr::pull(utr5end) + nnt_gene
} 


# Creates a matrix with number of columns that start at 5' UTR and finish at 50'th nt of the CDS
GetGeneDatamatrix5UTR <- function(gene, dataset, hdf5file, utr5_data, nnt_gene) {
  data_mat_all <- GetGeneDatamatrix(gene, dataset, hdf5file)
  #I'm changing n_left5 to just 1st nt of 5'UTR
  n_left5 <- Get5UTRstart(gene, utr5_data) # column to start from (5'end)
  n_right3 <- UTR5_length(gene, utr5_data) + nnt_gene # column to end with (3'end) 
  data_mat_5start <- data_mat_all[, n_left5 : n_right3]
  posn5start <- UTR5_5start(gene, utr5_data) # get utr5 start position (e.g. -250)
  posn5end <- UTR5_5end(gene, utr5_data) # get utr5 end position (e.g. -1) + nnt_gene
  data_mat_5start <- tibble::as_tibble(data_mat_5start, .name_repair="minimal")
  names(data_mat_5start) <- as.character(seq(from=posn5start, to=posn5end)) # add position as character 'name's to columns
  data_mat_5start <-tibble::add_column(data_mat_5start, read_length=10:50, .before=1) # minreadlength:maxreadlength
  return(data_mat_5start)
  #return(data_mat_5start, posn5start, posn5end)
}

# for one gene! (UTR5 + nnt_gene of CDS, with positions relative to CDS starting at 0, read_length is additional column)
utr5Matrix_1gene <- GetGeneDatamatrix5UTR(test_orfs[1], dataset, hdf5file, utr5_data, nnt_gene)
# > utr5Matrix_1gene
# # A tibble: 41 x 301
# read_length `-250` `-249` `-248` `-247` `-246` `-245` `-244` `-243` `-242`
# <int>  <int>  <int>  <int>  <int>  <int>  <int>  <int>  <int>  <int>
#   1        10      0      0      0      0      0      0      0      0      0
# 2          11      0      0      0      0      0      0      0      0      0
# 3          12      0      0      0      0      0      0      0      0      0

# utr5Matrix_allGene <- #
  # :s
  
# want to replicate this, but not with lapply.
# interesting <-lapply(test_orfs,
#                      function(gene)
#                        GetGeneDatamatrix5start(gene,
#                                                dataset,
#                                                hdf5file,
#                                                Get_Get5UTRstart(gene, gff_df),
#                                                UTR5full = UTR5_length(gene, gff_df),
#                                                nnt_gene)
# ) 


<<<<<<< HEAD
GetGeneDatamatrix5start <- function(gene, dataset, hdf5file, 
                                    Get5UTRstart, UTR5full, nnt_gene) {
  data_mat_all <- GetGeneDatamatrix(gene, dataset, hdf5file)
  #I'm changing the n_left to just 1st nt of 5'UTR
      n_left5 <- Get5UTRstart # column to start from (5'end)
    zeropad5_mat <- matrix(0, nrow = nrow(data_mat_all), ncol = 0)
   
  n_right3 <- UTR5full + nnt_gene - 1 # column to end with (3'end) 
  data_mat_5start <- data_mat_all[, n_left5:n_right3]
  return(cbind(zeropad5_mat, data_mat_5start))
}

interesting <-
  lapply(test_orfs[1],
         function(gene) 
           GetGeneDatamatrix5start(gene,
                                   dataset,
                                   hdf5file,
                                   Get5UTRstart(gene, gffdf = gff_df),
                                   UTR5full = UTR5_length(gene, gffdf = gff_df), 
                                   nnt_gene = 50)
  ) %>%
  Reduce("+", .) %>% # sums the list of data matrices
  TidyDatamatrix(startpos = -250 + 1, startlen = 10)


#maybe I can assume that the reduce bit is fine
#and then the TidyDataMatrix

#why is the startpos and startlen = 1?
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

### the functio
=======
# ------------------------------------------------------------------------------

##### BEGINNING OF WORKING SPACE

interesting <-lapply(test_orfs[1],
                     function(gene)
                       GetGeneDatamatrix5start(gene,
                                               dataset,
                                               hdf5file,
                                               Get_Get5UTRstart(gene, gff_df),
                                               UTR5full = UTR5_length(gene, gff_df),
                                               nnt_gene)
) 

###Instead of TidyDataMatrix: 
>>>>>>> 1aaa8cea6ea77b8a48e506e09048f49e5ac3be60

#all_genes_all_info <- as_tibble(interesting ) #I don't really get why that doesnt work, the only difference is that I am not specifying which gene even though it only contains 1 matrix. Am I opening up the matrix somehow with the second function? okay i think i know,the gene must be selected or otherwise it will turn all the lists into one tibble which isnt what we want 

#here i need to look into a function that separates separate matrix for each gene? or that at least treates each matrix separately 

try<-wtf %>% 
  set_colnames(positions) %>%
  gather(key="Position", value = "Counts") %>%
  group_by(Position) %>%
  summarise(Counts=sum(Counts))

try<-arrange(try,try$Position) ##ok slightly breaking here


##### END OF WORKING SPACE

interesting_try<-ggplot(
  data=try,) +
  geom_histogram(aes(x=Position, y=Counts), stat="identity") + 
  theme_light() +
  labs(y= "Read count", x = "Position (mRNA)")

ggsave("Another_one", device = "png")
#ok it's not the best plot ever but at least i know now what i need to do


### Potential functions?

Find() #returns the first element which matches the predicate (or the last element if right = TRUE). #That would be perfect if I wanted to for instnace only choose the first function at which it recognizes read counts


Position() #returns the position of the first element that matches the predicate (or the last element if right = TRUE).


integrate() #finds the area under the curve defined by f()
uniroot() #finds where f() hits zero
optimise() #finds the location of lowest (or highest) value of f()


<<<<<<< HEAD
### Playing
gene <- as.list(test_orfs[1])
gene %>%
GetGeneDatamatrix5start(gene,
                        dataset,
                        hdf5file,
                        Get5UTRstart= Get5UTRstart(gene, gffdf = gff_df),
                        UTR5full = UTR5_length(gene, gffdf = gff_df), 
                        nnt_gene = 50)
traceback()




 %>%
  Reduce("+", .) %>% # sums the list of data matrices
  TidyDatamatrix(startpos = -250 + 1, startlen = 10)

startpos:(startpos + ncol(data_mat) - 1)

play
=======
# 
# TidyDatamatrix <- function(data_mat, startpos = 1, startlen = 1) {
#   # CHECK startpos/off-by-one
#   positions <- startpos:(startpos + ncol(data_mat) - 1)
#   readlengths <- startlen:(startlen + nrow(data_mat) - 1)
#   data_mat %>%
#     set_colnames(positions) %>%
#     as_tibble() %>%
#     mutate(ReadLen = readlengths) %>%
#     gather(-ReadLen, key = "Pos", value = "Counts", convert = FALSE) %>%
#     mutate(Pos = as.integer(Pos), Counts = as.integer(Counts))
# }
>>>>>>> 1aaa8cea6ea77b8a48e506e09048f49e5ac3be60


##Selecting 5'UTR start position required for plotting
#the function would identify the first position which identifies read counts and from that (or    this position+ for instance 10 nt just for visualization)

