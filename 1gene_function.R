#### Starting with one gene: YCR012W (PGK1)

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

###Set up:

# prepare files, opens hdf5 file connection
dataset <- "G-Sc_2014"
hd_file <- "G-Sc_2014/output/WTnone/WTnone.h5"
hdf5file <- rhdf5::H5Fopen(hd_file) # filehandle for the h5 file

#Initial set ups
gene1 <- "YCR012W"
test_orfs <- c("YCR012W","YEL009C","YOR303W","YOL130W","YGR094W")
nnt_gene<- 50
startpos <-250
startlen <- 10
orf_gff_file <- "G-Sc_2014/input/yeast_CDS_w_250utrs.gff3"
gff_df <- readGFFAsDf(orf_gff_file)

##Functions

#This will be used to create a matrix from hdf5file
GetGeneDatamatrix <- function(gene, dataset, hdf5file) {
  hdf5file %>%
    rhdf5::H5Dopen(
      name = paste0("/", gene, "/", dataset, "/reads/data")
    ) %>%
    rhdf5::H5Dread() %>%
    return()
}

#Function to create a gff table  
readGFFAsDf <- purrr::compose(
  rtracklayer::readGFFAsGRanges,
  data.frame, 
  as_tibble,
  .dir = "forward" # functions called from left to right
)

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
#Currently replaced by TidyDataMatrix 
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

#plotting function
plot<- function(x, gene) {
  my_title <-title(gene)  #it wooooooorks!
  text_AUG <- textGrob("AUG", gp=gpar(fontsize=13, fontface="bold")) #to place text annotation
  #the actual plotting 
  plotted_UTR <- ggplot(x) +
    geom_density(aes(x=Pos, y=Counts), stat="identity") +
    scale_x_continuous(limits = c(-250,50), expand = c(0, 0)) +
    scale_y_continuous(expand = c(0,0)) +
    labs(y= "Read count", x = "Position") +
    ggtitle(my_title) + 
    annotate("text", x = 0, y = -3, label = "Some text") +
    # coord_cartesian(clip = "off") + 
    # annotation_custom(text_AUG,xmin=0,xmax=0,ymin=0,ymax=5) + 
    theme_classic()
}

##Tried and doesnt work:
 #annotate("text", x = 0, y = -3, label = "Some text")
  #ok i know why the annotation custom doesnt work-that's because y axis is different for each graph

title <- function(gene) {
  paste("Ribosome footprint density of ", gene , sep= "")
}

#final function from raw data processing to data visualization
UTR5_plot <- function(gene) {
  lapply(gene,
         function(gene) 
           GetGeneDatamatrix5UTR(gene,
                                 dataset,
                                 hdf5file,
                                 x=gff_df,
                                 nnt_gene = nnt_gene)
  )%>%
    Reduce("+", .) %>% # sums the list of data matrices
    TidyDatamatrix(startpos = -250, startlen = 10) %>%
  plot(. , gene) %>%
  return()
}

#the final function i'd use :)
final_function <- function(x){
 purrr::map(x, UTR5_plot)
}

##working space:

UTR5_plot_test <- function(gene) {
  lapply(gene,
         function(gene) 
           GetGeneDatamatrix5UTR(gene,
                                 dataset,
                                 hdf5file,
                                 x=gff_df,
                                 nnt_gene = nnt_gene)
  )%>%
    Reduce("+", .) %>% # sums the list of data matrices
    TidyDatamatrix(startpos = -250, startlen = 10)
}
xxx<-UTR5_plot_test(test_orfs[1]) #no AUGs for any of them except 1 if done that way 
plotx<-plot(xxx,test_orfs[1])
plotx
