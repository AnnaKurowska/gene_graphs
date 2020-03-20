#### Starting with one gene: YCR012W (PGK1)

############################################################################################
                                        ##set up##
setwd("/Users/Ania/Desktop/Szkoła/4th year/Dissertation/gene_graphs/")

##Libraries 

# need this
library(rhdf5)
library(tidyr)
library(dplyr)
library(magrittr)
library(purrr)
library(ggplot2)
library(ggpubr)
library(tibble)
library(grid)
 #for the threshold function

# library(expss) #for count_col_if
# library(HH) #  for position() I probably won't need any of them 

###Set up:
#WT_none <- "G-Sc_2014/output/WTnone/WTnone.h5"
WT_3AT <- "G-Sc_2014/output/WT3AT/WT3AT.h5"
#WT_CHX <- "G-Sc_2014/output/WTCHX/WTCHX.h5"

# prepare files, opens hdf5 file connection
dataset <- "G-Sc_2014"
hd_file <- WT_3AT
hdf5file <- rhdf5::H5Fopen(hd_file) # filehandle for the h5 file

#Initial set ups
test_orfs <- c("YCR012W","YEL009C","YOR303W","YOL130W","YGR094W")
nnt_gene<- 50
startpos <-250
startlen <- 10
temporary_length <- 20 #for zooming in bit
min_read_count <- 10
orf_gff_file <- "G-Sc_2014/input/yeast_CDS_w_250utrs.gff3"

#Function to create a gff table  
readGFFAsDf <- purrr::compose(
  rtracklayer::readGFFAsGRanges,
  data.frame, 
  as_tibble,
  .dir = "forward" # functions called from left to right
)

gff_df <- readGFFAsDf(orf_gff_file) # need to run this first readGFFAsDf (below)
############################################################################################
                                    ##Functions##
#This will be used to create a matrix from hdf5file
GetGeneDatamatrix <- function(gene, dataset, hdf5file) {
  hdf5file %>%
    rhdf5::H5Dopen(
      name = paste0("/", gene, "/", dataset, "/reads/data")
    ) %>%
    rhdf5::H5Dread() %>%
    return()
}

# Takes value for start position of 5'UTR for matrix creation (no negative values e.g.: 1)
Get5UTRstart <- function(gene, x) {
  x %>% 
    dplyr::filter(Name == gene, type=="UTR5" ) %>% 
    dplyr::pull(start)
}

# Takes value for 3' end of the table (includes whole 5'UTR +50 nt of CDS)
CDS3_end <- function(gene, x) {
  x %>% 
    dplyr::filter(Name == gene, type=="UTR5") %>% 
    dplyr::pull(end) + nnt_gene
} 


# Takes value for length of 5'UTR  
UTR5_length <- function(gene, x) {
  x %>% 
    dplyr::filter(Name == gene, type=="UTR5" ) %>% 
    dplyr::pull(width)
} 

# Creates a matrix with number of columns that start at 5' UTR and finish at 50'th nt of the CDS
GetGeneDatamatrix5UTR <- function(gene, dataset, hdf5file, x, nnt_gene) {
  data_mat_all <- GetGeneDatamatrix(gene, dataset, hdf5file)
  n_left5 <- Get5UTRstart(gene, x) # column to start from (5'end)
  n_right3 <- UTR5_length(gene, x) + nnt_gene # column to end with (3'end) 
  data_mat_5start <- data_mat_all[, n_left5 : n_right3]
  data_mat_5start <- tibble::as_tibble(data_mat_5start, .name_repair="minimal")
  return(data_mat_5start)
  #return(data_mat_5start, posn5start, posn5end)
}

# GetPosnCountOutput <- function(x){
#   output_thing <- x %>% 
#     gather(-read_length, key="Position", value = "Counts") %>%
#     group_by(Position) %>%
#     summarise(Counts=sum(Counts)) %>%
#     arrange(as.integer(Position))
#   return(output_thing)
# }

TidyDatamatrix <- function(x,startpos = 1, startlen = 1,gene) {
  # CHECK startpos/off-by-one
  positions <- startpos:(startpos + ncol(x) - 1)
  readlengths <- startlen:(startlen + nrow(x) - 1)
  x %>%
    set_colnames(positions) %>%
    as_tibble() %>% 
    # names(x) <- paste(gene) %>%  how to do that?
    mutate(ReadLen = readlengths) %>%
    gather(-ReadLen, key = "Pos", value = "Counts", convert = FALSE) %>%
    mutate(Pos = as.integer(Pos), Counts = as.integer(Counts)) %>%
    group_by(Pos) %>%
    summarise(Counts=sum(Counts))
}

#final function from raw data processing to data visualization
UTR5_table <- function(gene) {
  lapply(gene,
         function(gene) 
           GetGeneDatamatrix5UTR(gene,
                                 dataset,
                                 hdf5file,
                                 x=gff_df,
                                 nnt_gene = nnt_gene)
  )%>% #names() <-paste0(gene, SAMPLETBC)
    Reduce("+", .) %>% # sums the list of data matrices
    TidyDatamatrix(startpos = -250, startlen = 10) %>%
    # uAUG_efficiency %>%  optional to put it here
  return()
}


#the final function i'd use :) iteration of UTR5_table 
final_function_table <- function(x){
 purrr::map(x, UTR5_table)  ### as_tibble(.name_repair = "unique)
}

output_orfs<-final_function_table(test_orfs) #it works just fine 


############################################################################################                                               ##plotting##

#main plotting function 1gene-1 plot (iterated later on in plotting_multiple)
plotting_5UTR<- function(input_data) {
  text_AUG <- textGrob("AUG", gp=gpar(fontsize=13, fontface="bold")) #to place text annotation

  #the actual plotting
  plotted_UTR <- ggplot(input_data) +
    geom_density(aes(x=Pos, y=Counts), stat="identity") +
    scale_x_continuous(limits = c(-250,50), expand = c(0, 0)) +
    scale_y_continuous(expand = c(0,0)) +
    labs(y= "Read count", x = "Position") +
    ggtitle(paste("Ribosome footprint density of ", names(input_data) , sep= "")) +
    # coord_cartesian(clip = "off") +
    # annotation_custom(text_AUG,xmin=0,xmax=0,ymin=0,ymax=5) +
    theme_classic()
}
                            ##multiple plots##

plotting_multiple <- function(input_data) {
  purrr :: map(input_data, plotting_5UTR) %>%
    ggarrange(plotlist = .) %>%  ##, common.legend = (maybe to establish a name)
    return()
    #ggsave(file.path(paste0("Plot of ", gene)), device = "jpg")
}

####################################################################################
                                    ##2 in one##

#
# hd_file <- WT_3AT
# hdf5file <- rhdf5::H5Fopen(hd_file) # filehandle for the h5 file
#
# output_orfs_3AT <-final_function_table(test_orfs)
#
# ######

# two_in_one_plots(output_orfs[1], output_orfs_3AT[1]) # they're still lists and we want tibbles, otherwise it doesnt work
#
# two_in_one_plots <- function(sample1, sample2) {
#   more_plots <- ggplot(sample1) +
#     geom_density(aes(x=Pos, y=Counts), stat="identity") +
#     geom_density(data = sample2, aes(x=Pos, y= Counts), stat="identity", color = "red") +
#     scale_x_continuous(limits = c(-100,50), expand = c(0, 0)) +
#     scale_y_continuous(expand = c(0,0)) +
#     labs(y= "Read count", x = "Position") +
#     theme_classic()
# }

###################################################################################
                            #A site displacement slot#

###
CalcAsiteFixedOneLength <- function(reads_pos_length, min_read_length,
                                    read_length, asite_disp) {
  # Calculate read A-site using a fixed displacement for a single read length
  length_row_choose <- read_length - min_read_length + 1
  reads_pos_length[length_row_choose, ] %>%
    dplyr::lag(n = asite_disp, default = 0)
}

something <- CalcAsiteFixedOneLength(reads_pos_length = GetGeneDatamatrix(gene = test_orfs[1],dataset = dataset,hdf5file = hdf5file), min_read_length = 10, read_length = 28, asite_disp = 15)
# > str(something)
# num [1:1751] 0 0 0 0 0 0 0 0 0 0 ...

# > something
# [1]    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
# [20]    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
# [39]    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
# [58]    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
# [77]    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
# [96]    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
# [115]    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
# [134]    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
# [153]    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
# [172]    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    1    0    0
# [191]    0    0    0    0    0    0    1    0    0    1    0    0    2    0    0    0    3    0    0
# [210]    1    0    9    1    0    1    0    0    4    5    0    0    0    0    0    0    0    0    0
# [229]    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
# [248]    0    0    0    0    0    0   80    8    0    6    2    0    6    1    0    9    2    0   59
# [267]   18    0    6    6    3  129    6    1  278   20    0   40   17    3   42    3    0    1    1
# [286]    1   19    0    0    0    4    9   36   14    0   42    6   33   25    2    0   12    1    0
# [305]   22    1    0    1    0    0  162    0    3   71   51    0   73   53   18   55   25    1   37
# [324]    3    0  502   65   10   83   54    6   10   20    7   19    2    0   53    1    0   11    0
# [343]    1   27   53    1  279   30    1  283   13    0    0   23    0    0    3    0   10    0    1
# [362]    7    3    0   62   17    2  120   62    4  107   21   21   10    5    2   40  127    0  101
# [381]   12    5   91   53    1   33   34    2    8   11   12  139   28   16   13    9    0   87   32
# [400]    0  114   18    2  352   30   13 2532   28   13   86  137    2  273   70   38   61   23    1

###
CalcAsiteFixed <- function(reads_pos_length, min_read_length,
                           asite_disp_length = data.frame(
                             read_length = c(28, 29, 30),
                             asite_disp = c(15, 15, 15)
                           ),
                           colsum_out = TRUE) {
  # Calculate read A-site using a fixed displacement for fixed read lengths
  npos <- ncol(reads_pos_length)
  Asite_counts_bylength <-
    purrr::map2(
      asite_disp_length$read_length, asite_disp_length$asite_disp,
      function(read_length, asite_disp) {
        CalcAsiteFixedOneLength(
          reads_pos_length,
          min_read_length,
          read_length,
          asite_disp
        )
      }
    )
  if (colsum_out) {
    Asite_counts <- purrr::reduce(Asite_counts_bylength, `+`)
    return(Asite_counts)
  } else {
    # this has only as many columns as asite_disp_length,
    # probably LESS than data_mat
    Asite_counts_bylengthmat <- unlist(Asite_counts_bylength) %>%
      matrix(ncol = npos, byrow = TRUE)
    return(Asite_counts_bylengthmat)
  }
}

CalcAsiteFixed(reads_pos_length =GetGeneDatamatrix(gene = test_orfs[1],dataset = dataset,hdf5file = hdf5file), min_read_length = 10 , asite_disp_length =data.frame(
  read_length = c(28, 29, 30),
  asite_disp = c(15, 15, 15)
) , colsum_out = TRUE)

# [1]    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
# [17]    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
# [33]    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
# [49]    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
# [65]    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
#.....
# [ reached getOption("max.print") -- omitted 751 entries ]

###
SnapToCodon <- function(x, left, right, snapdisp=0L) {
  # snap nucleotide-aligned reads to codon position
  #   x:     vector
  #   left:  integer for starting position, frame 0
  #   right: integer for ending position
  #   snapdisp: integer any additional displacement in the snapping
  RcppRoll::roll_suml(x[(left:right) + snapdisp], n=1L, by=1L, fill = NULL)
}  ##I changed n= and by= from 3L to 1L bc I want NT resolution 

left = Get5UTRstart(test_orfs[1], gff_df) #finds left and right value for each gene 
right = CDS3_end(test_orfs[1], gff_df)

# Get5UTRstart <- function(gene, x) {
#   x %>% 
#     dplyr::filter(Name == gene, type=="UTR5" ) %>% 
#     dplyr::pull(start)
# }
# 
# # Takes value for 3' end of the table (includes whole 5'UTR +50 nt of CDS)
# CDS3_end <- function(gene, x) {
#   x %>% 
#     dplyr::filter(Name == gene, type=="UTR5") %>% 
#     dplyr::pull(end) + nnt_gene
# } 

###
GetGeneCodonPosReads1dsnap <- function(gene, dataset, hdf5file, left, right, 
                                       min_read_length, 
                                       asite_disp_length = data.frame(
                                         read_length = c(28, 29, 30),
                                         asite_disp = c(15, 15, 15)
                                       ), 
                                       snapdisp=0L) {
  reads_pos_length <- GetGeneDatamatrix(gene, dataset, hdf5file) # Get the matrix of read counts
  reads_asitepos <- CalcAsiteFixed(
    reads_pos_length, min_read_length,
    asite_disp_length
  )
  SnapToCodon(reads_asitepos,left,right,snapdisp)
}

final_A_mapped <- GetGeneCodonPosReads1dsnap(gene = test_orfs[1], dataset = dataset ,hdf5file = hdf5file ,left =left , right = right , min_read_length = 10 , asite_disp_length = data.frame(
  read_length = c(28, 29, 30),
  asite_disp = c(15, 15, 15)), snapdisp = 0L)

> final_A_mapped
# [1]   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
# [22]   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0
# [43]   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1
# [64]   1   0   1   1  10   4   9  15   3  89   3   0   0   0   0   0   0   0   1   0   0
# [85] 142  20  19  13  82  15 225 409  79  51   7  28  12 174  85 225
# > 
# For RPF datasets, generate codon-based position-specific reads

    # Get codon-based position-specific reads for each gene, in a tibble
    # reads_per_codon_etc <- tibble(gene=test_orfs[1]) %>%
Counts_Asite_mapped  <- map(test_orfs, ~GetGeneCodonPosReads1dsnap(
  .,
  dataset,
  hdf5file,
  left = left,
  right = right,
  min_read_length = 10,
  asite_disp_length = data.frame(
    read_length = c(28, 29, 30),
    asite_disp = c(15, 15, 15)
  ))) 

%>%  
  as_tibble(.name_repair = "unique")
 
  function(x,startpos = 1, startlen = 1,gene) {
    # CHECK startpos/off-by-one
    positions <- startpos:(startpos + ncol(x) - 1)
    readlengths <- startlen:(startlen + nrow(x) - 1)
    x %>%
      set_colnames(positions) %>%
      as_tibble() %>% 
      # names(x) <- paste(gene) %>%  how to do that?
      mutate(ReadLen = readlengths) %>%
      gather(-ReadLen, key = "Pos", value = "Counts", convert = FALSE) %>%
      mutate(Pos = as.integer(Pos), Counts = as.integer(Counts)) %>%
      group_by(Pos) %>%
      summarise(Counts=sum(Counts))
  }
   # as_tibble(.name_repair = "unique") %>%
  # set_colnames(paste(test_orfs))  %>%
  # cbind(seq(from = -250, to = 49, by = 1))




# tibble1 <- output_orfs[[1]]
# A_site_and_counts <- cbind(Counts_Asite_mapped, tibble1)
# plotting_5Asite<- function(input_data) {
#   text_AUG <- textGrob("AUG", gp=gpar(fontsize=13, fontface="bold")) #to place text annotation
#   
#   #the actual plotting
#   plotted_UTR <- ggplot(input_data) +
#     geom_density(aes(x=Pos, y=CountsAsite), stat="identity") +
#     scale_x_continuous(limits = c(-250,50), expand = c(0, 0)) +
#     scale_y_continuous(expand = c(0,0)) +
#     labs(y= "Read count", x = "Position") +
#     ggtitle(paste("Ribosome footprint density of ", names(input_data) , sep= "")) +
#     # coord_cartesian(clip = "off") +
#     # annotation_custom(text_AUG,xmin=0,xmax=0,ymin=0,ymax=5) +
#     theme_classic() %>%
#     return()
# }
# A_site_and_counts_plot <-plotting_5Asite(A_site_and_counts)
# comparison_plot <-plotting_5UTR(tibble1)
# ggarrange(A_site_and_counts_plot,comparison_plot)

####################################################################################
                                      ##Zooming in ## 

# test object with 1 gene
#tibble1 <- output_orfs[[1]]
# > str(tibble1)
# Classes ‘tbl_df’, ‘tbl’ and 'data.frame':	300 obs. of  2 variables:
#   $ Pos   : int  -250 -249 -248 -247 -246 -245 -244 -243 -242 -241 ...
# $ Counts: int  0 0 0 0 0 0 0 0 0 0 ...

# tibble5 <- output_orfs[[5]] i've created more tables to play with 


##### Function for the upstream codons 

# finding_uAUG_beginning <- function(datatibble) {
#   datatibble %>% 
#     dplyr :: filter(Counts > 0) %>%
#     min() %>%
#     return()
# }

###i was playing with tibble5 (instead of tbl1) bc it has much lower count frequency, i encountered the following issue: 
#for finding_uAUG_beginning(tibble5) if filter(Counts > 4) bc the max count value of tbl5 is 4 so that's a good question how i can tackle that. 
# Error in FUN(X[[i]], ...) : 
#   only defined on a data frame with all numeric variables 

# to run: 
#finding_uAUG_beginning(tibble1)
# example run, set to uAUG_beginning:
#uAUG_beginning <- finding_uAUG_beginning(tibble1)

#finding_uAUG_ending <- function(uAUG_start) {
#   uAUG_start + temporary_length %>%
#     return()
# }

# to run: 
#finding_uAUG_ending(uAUG_beginning)
# example run, set to uAUG_ending:
#uAUG_ending <- finding_uAUG_ending(uAUG_beginning)

# sum_uAUG <- function(datatibble, uAUG_start, uAUG_end){
#   datatibble %>%
#     filter(Pos >= uAUG_start, Pos <= uAUG_end) %>%
#     summarise(sum_read_counts_uAUG=sum(Counts)) %>%
#     return()
# }
# 
# # to run:
# #sum_uAUG(tibble1, uAUG_beginning, uAUG_ending)
# # example run set to sum_uAUG_result
# #sum_uAUG_result <- sum_uAUG(tibble1, uAUG_beginning, uAUG_ending)
# 
# 
# ####AUG
# finding_AUG_beginning <- function(datatibble) {
#   datatibble %>%
#     dplyr::select(Pos) %>%
#     filter(Pos == -10) %>%  #I changed the value to 6 so altogether now it's 12nt region
#     pull
# }
# 
# # to run:
# #finding_AUG_beginning(tibble1)
# # example run, set to uAUG_beginning:
# #AUG_beginning <- finding_AUG_beginning(tibble1)
# 
# finding_AUG_ending <- function(AUG_start) {
#     AUG_start + temporary_length %>%
#     return()
# }
# 
# # to run:
# #finding_AUG_ending(AUG_beginning)
# # example run, set to uAUG_ending:
# #AUG_ending <- finding_AUG_ending(AUG_beginning)
# 
# sum_AUG <- function(datatibble, AUG_start, AUG_end){
#   datatibble %>%
#     filter(Pos >= AUG_start, Pos <= AUG_end) %>%
#     summarise(sum_read_counts_AUG=sum(Counts)) %>%
#     return()
# }
# # to run:
# #sum_AUG(tibble1, AUG_beginning, AUG_ending )
# # example run, set to uAUG_ending:
# #sum_AUG_result <-sum_AUG(tibble1, AUG_beginning, AUG_ending )
# 
# ##efficiency
# upstream_efficiency <- function(uAUG_sum_result, AUG_sum_result) {
#   (uAUG_sum_result/AUG_sum_result *100) %>%
#     return()
# }
# 
# #upstream_efficiency(sum_uAUG_result/sum_AUG_result)
# 
# ##final function
# 
# uAUG_efficiency <- function(datatibble, uAUG_start,uAUG_end) {
# #finding uAUG and AUG endpoints
#   AUG_start <- -10
#   AUG_end <- AUG_start + temporary_length
#   # AUG_start <-finding_AUG_beginning(datatibble)
#   # AUG_end <- finding_AUG_ending(AUG_start)
#   # uAUG_start <-finding_uAUG_beginning(datatibble)
#   # uAUG_end <- finding_uAUG_ending(uAUG_start)
# 
#   #summing uAUG and AUG region
#   final_sum_uAUG <- sum_uAUG(datatibble, uAUG_start, uAUG_end )
#   final_sum_AUG <- sum_AUG(datatibble, AUG_start, AUG_end )
# 
#   #calculating the efficiency
#   final_result_of_efficiency <- upstream_efficiency(final_sum_uAUG, final_sum_AUG)
#   paste0(final_result_of_efficiency)
# }
# 
# ##
# uAUG_efficiency_many_genes <- function(datatibble, uAUG_start, uAUG_end) {
#   map(datatibble, uAUG_efficiency, uAUG_start, uAUG_end)
# }
# ## example:
# saved <- uAUG_efficiency_many_genes(output_orfs,uAUG_start = -Inf,uAUG_end = -11 ) %>%
#   as_tibble(.name_repair = "minimal")

###########
### Creating a table for multiple positions for multiple genes

#0:-15 nt
#0:-60 nt
#60:-120 nt
#0:-250 nt

#creating adjusted regions to sum read counts 
sum_uAUG <- function(datatibble, uAUG_start, uAUG_end){
  datatibble %>%
    filter(Pos >= uAUG_start, Pos <= uAUG_end) %>%
    summarise(sum_read_counts_uAUG=sum(Counts)) %>%
    pull
}

#for multiple genes 
sum_uAUG_multiple<- function(genes, gene_names, uAUG_start, uAUG_end) {
  #function 
  map(genes, sum_uAUG, uAUG_start, uAUG_end) %>%
    tibble::enframe(name = NULL)  %>%
    set_rownames(., paste0(gene_names)) %>%
  #names of columns and rows
  return()
}

all_regions_together <- function(genes, gene_names) {
  region1 <-sum_uAUG_multiple(genes, gene_names,  uAUG_start = -15, uAUG_end = 0) 
  region2 <-sum_uAUG_multiple(genes, gene_names,  uAUG_start = -60, uAUG_end = 0)
  region3 <-sum_uAUG_multiple(genes, gene_names,  uAUG_start = -120, uAUG_end = -60)
  region4 <-sum_uAUG_multiple(genes, gene_names, uAUG_start = -250, uAUG_end = 0)
  
  cbind(region1, region2, region3, region4) %>%
    set_colnames(c("-15:0", "-60:0", "-120:-60", "-250:0")) %>%
    return()
}

all_regions_together(genes = output_orfs,gene_names = test_orfs)

################################ seems we're not using it at the end but still useful 
sliding_windows <- function(tibble) {
  RcppRoll::roll_sum(tibble$Counts, n =30, by = 30) %>% 
    tibble::enframe(name = NULL) %>%
     t()
}

sliding_windows_multiple <- function(tibble, genes){
  names(tibble) <- genes
  joined <-map(tibble, sliding_windows ) 
  col_names <- seq(from = -250, to = 20, by = 30)
  
  bind_rows(joined) %>%
    set_rownames( value = c(seq(from = -250, to = 20, by = 30)))  %>%
    t() %>%
    return()  ###add a column instead of the set_rownames and you want these to be ranges (from -250 to -230)
}

xxxx <-sliding_windows_multiple(output_orfs, test_orfs)
##########################################################################################
                                  ##AUG annotation issue 

# ##working space: for the AUG annotation issue 
# 
# UTR5_plot_test <- function(gene) {
#   lapply(gene,
#          function(gene) 
#            GetGeneDatamatrix5UTR(gene,
#                                  dataset,
#                                  hdf5file,
#                                  x=gff_df,
#                                  nnt_gene = nnt_gene)
#   )%>%
#     Reduce("+", .) %>% # sums the list of data matrices
#     TidyDatamatrix(startpos = -250, startlen = 10)
# }
# xxx<-UTR5_plot_test(test_orfs[1]) #no AUGs for any of them except 1 if done that way 
# plotx<-plot(xxx,test_orfs[1])
# plotx
