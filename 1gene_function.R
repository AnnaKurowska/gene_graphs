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
WT_none <- "G-Sc_2014/output/WTnone/WTnone.h5"
WT_3AT <- "G-Sc_2014/output/WT3AT/WT3AT.h5"
WT_CHX <- "G-Sc_2014/output/WTCHX/WTCHX.h5"

# prepare files, opens hdf5 file connection
dataset <- "G-Sc_2014"
hd_file <- WT_3AT
hdf5file <- rhdf5::H5Fopen(hd_file) # filehandle for the h5 file

#Initial set ups
test_orfs <- c("YCR012W","YEL009C","YOR303W","YOL130W","YGR094W")
nnt_gene<- 50
startpos <-250
startlen <- 10
temporary_length <- 12 #for zooming in bit
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
  return()
}


#the final function i'd use :) iteration of UTR5_table 
final_function_table <- function(x){
 purrr::map(x, UTR5_table)  ### as_tibble(.name_repair = "unique)
}

output_orfs<-final_function_table(test_orfs) #it works just fine 


############################################################################################                                        ##plotting##

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

###ok so that's the final function for plotting all graphs from a list of tibbles (output from )
    #still dont know how to name each graph and give a title to the plot overall
plotting_multiple <- function(input_data) {
  purrr :: map(input_data, plotting_5UTR) %>%
    ggarrange(plotlist = .) %>%  ##, common.legend = (maybe to establish a name)
    return()
    #ggsave(file.path(paste0("Plot of ", gene)), device = "jpg")
}

xxx<-plotting_multiple(output_orfs) 
#issue with the name now..same for all, i think it means we need to name them earlier on so we know what's what. Why does it only use the YCR012W?

ggsave("multiple", device = "jpg")

# ###### 2 in one
# 
# # for WT_3AT
# 
# hd_file <- WT_3AT
# hdf5file <- rhdf5::H5Fopen(hd_file) # filehandle for the h5 file
# 
# output_orfs_3AT <-final_function_table(test_orfs)
# 
# ######
# #so now we need to create a plot including both 
# ##we need to select each tibble only (can use map again)
# 
# 
# ###okay i absolutely can't be bothered to work on it anymore, im gonna move to the zooming in bit
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

##########################################################################################
                                      ##Zooming in ## 

# test object with 1 gene
tibble1 <- output_orfs[[1]]
# > str(tibble1)
# Classes ‘tbl_df’, ‘tbl’ and 'data.frame':	300 obs. of  2 variables:
#   $ Pos   : int  -250 -249 -248 -247 -246 -245 -244 -243 -242 -241 ...
# $ Counts: int  0 0 0 0 0 0 0 0 0 0 ...
tibble5 <- output_orfs[[5]]

#test run:
uAUG_efficiency(tibble5)
# example run, for the final function:
efficiency_final_full_result <- uAUG_efficiency(tibble1)
"#The efficiency of upstream translation initiation is  0.119617224880383 %"

#test run for multiple genes:
map(output_orfs, uAUG_efficiency)
# [[1]]
# [1] "The efficiency of upstream translation initiation is  0.119617224880383 %"
# 
# [[2]]
# [1] "The efficiency of upstream translation initiation is  Inf %"
# 
# [[3]]
# [1] "The efficiency of upstream translation initiation is  266.666666666667 %"
# 
# [[4]]
# [1] "The efficiency of upstream translation initiation is  113.333333333333 %"
# 
# [[5]]
# [1] "The efficiency of upstream translation initiation is  Inf %"
###needs thresholds, [1] encountered one single read count and it started counting from it
###needs a solution for Inf as well (I think it's division by 0)

uAUG_efficiency <- function(datatibble) {
  #finding uAUG and AUG start and end sites for counting
  uAUG_start <-finding_uAUG_beginning(datatibble)
  uAUG_end <- finding_uAUG_ending(uAUG_start)
  AUG_start <-finding_AUG_beginning(datatibble)
  AUG_end <- finding_AUG_ending(AUG_start)

  #summing uAUG and AUG region
  final_sum_uAUG <- sum_uAUG(datatibble, uAUG_start, uAUG_end )
  final_sum_AUG <- sum_AUG(datatibble, AUG_start, AUG_end )

  #calculating the efficiency
  upstream_efficiency(final_sum_uAUG, final_sum_AUG) %>%
    return()
}

##### Function for the upstream codons 
finding_uAUG_beginning <- function(datatibble) {
  datatibble %>% 
    dplyr :: filter(Counts > 0) %>%
    min() %>%
    return()
}

#for finding_uAUG_beginning(tibble5) if filter(Counts > 4) bc the max count value of tbl5 is 4 so that's a good question how i can tackle that. 
# Error in FUN(X[[i]], ...) : 
#   only defined on a data frame with all numeric variables 

finding_uAUG_beginning(tibble1) #works for a single value
map(output_orfs, finding_uAUG_beginning)# doesn't work if it's larger than 4
##deciphering:


# to run: 
finding_uAUG_beginning(tibble1)
# example run, set to uAUG_beginning:
uAUG_beginning <- finding_uAUG_beginning(tibble1)

finding_uAUG_ending <- function(uAUG_start) {-
  uAUG_start + temporary_length %>%
    return()
}
# to run: 
finding_uAUG_ending(uAUG_beginning)
# example run, set to uAUG_ending:
uAUG_ending <- finding_uAUG_ending(uAUG_beginning)

sum_uAUG <- function(datatibble, uAUG_start, uAUG_end){
  datatibble %>%
    filter(Pos >= uAUG_start, Pos <= uAUG_end) %>%
    summarise(sum_read_counts_uAUG=sum(Counts)) %>%
    return()
}
# to run: 
sum_uAUG(tibble1, uAUG_beginning, uAUG_ending)
# example run set to sum_uAUG_result
sum_uAUG_result <- sum_uAUG(tibble1, uAUG_beginning, uAUG_ending)

############################################################################################
                                 ##Working space for tibble5##

##### Function for the AUG codons 
finding_AUG_beginning <- function(datatibble) {
  datatibble %>% 
    dplyr::select(Pos) %>%
    filter(Pos == -6) %>%  #I changed the value to 6 so altogether now it's 12nt region
    pull 
}

#so i'm first calculating the efficiency of AUG initiation 
#start of the search, i could be just looking at a specific region flanking the AUG site, maybe it could be maybe 24 NTs?

# to run: 
finding_AUG_beginning(tibble1)
# example run, set to AUG_beginning:
AUG_beginning <- finding_AUG_beginning(tibble5)
## we don't want a tibble, it stops it from working 
# A tibble: 1 x 1
# Pos
# <int>
#   1     0

# AUG_start_tibble5 and AUG_end_tibble5 and AUG_sum_tibble5:
AUG_sum_tibble5 <- sum_AUG(tibble5, AUG_start_tibble5, AUG_end_tibble5)
# A tibble: 1 x 1
# sum_read_counts_AUG
# <int>
#   1                   

##idea: find different positions at which read count > 0, calculate the efficiency of that site relative to the AUG site, (which must be the same for AUG and each non-AUG), start counting from counts>0: if larger than some % threshold then calculate the efficiency, if smaller just skip, do that for each counts > 0 result, the output will be only done for sites with a larger threshold 

#step: find each counts > 0 and calculate the % for it 

#you're getting a table with only counts > 0 values,
finding_uAUG_beginning2 <- function(datatibble) {
  datatibble %>% 
    dplyr :: filter(Counts > 0) %>%
    pull(Pos)
}

#find beginning for each count > 0 site
AUG_beginning_tibble5 <- finding_uAUG_beginning2(tibble5)
# > AUG_beginning_tibble5
#[1] -31 -25 -22 -17 -16 -14 -13 -10  -8   0  24  27  32  41  42  49

#find end for each count > 0 site
AUG_end_tibble5<-finding_uAUG_ending(AUG_beginning_tibble5)
#[1] -19 -13 -10  -5  -4  -2  -1   2   4  12  36  39  44  53  54  61



#calculate the sum for each counts>
map(AUG_beginning_tibble5,sum_uAUG )


#then calculate the percentage relative to AUG 


############################################################################################


finding_AUG_ending <- function(AUG_start) {
  AUG_start + temporary_length %>% 
    return()
}

# to run: 
finding_AUG_ending(AUG_beginning)
# example run, set to AUG_ending:
AUG_ending <- finding_AUG_ending(AUG_beginning)

sum_AUG <- function(datatibble, AUG_start, AUG_end){
  datatibble %>%
    filter(Pos >= AUG_start, Pos <= AUG_end) %>%
  summarise(sum_read_counts_AUG=sum(Counts)) %>%
    return()
}
# to run: 
sum_AUG(tibble1, AUG_beginning, AUG_ending)
# example run, set to sum_AUG:
sum_AUG_result <- sum_AUG(tibble1, AUG_beginning, AUG_ending)
# # A tibble: 1 x 1
# sum_read_counts_AUG
# <int>
#   1                 836

##### uAUG efficiency relative to AUG codon 
upstream_efficiency <- function(sum_uAUG, sum_AUG) {
  efficiency <- ((sum_uAUG/sum_AUG) * 100) 
    paste0("The efficiency of upstream translation initiation is  ", efficiency, " %" ) %>%
      return()
} 
#to run
upstream_efficiency(sum_uAUG_result,sum_AUG_result)
# example run, set to upstream_efficiency:
efficiency_value <- upstream_efficiency(sum_uAUG_result,sum_AUG_result)

##
#threshold

###so this may actually be very different, i may not want to find a single value, but instead 1) use a threshold abd 2) 

                                      ##Zooming in ##
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
