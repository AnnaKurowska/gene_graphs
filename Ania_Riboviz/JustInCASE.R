
############################################################################################
##set up##
# setwd("/Users/Ania/Desktop/Szkoła/4th year/Dissertation/gene_graphs/")
library(here)
here()
# "/Users/Ania/Desktop/planets" is that the issue for Mrkd?
##Libraries 


# if (!require(rhdf5){ install.packages("rhdf5")
#   library(rhdf5)
# }
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
library(stringr)
library(ape)
library(ggrepel)

#set up for Guydosh
WT_none <- "G-Sc_2014/output/WTnone/WTnone.h5"
hd_file_none <- WT_none
hdf5file_none <- rhdf5::H5Fopen(hd_file_none) # filehandle for the h5 file

## WT_3AT
WT_3AT <- "G-Sc_2014/output/WT3AT/WT3AT.h5"
hd_file_3AT <- WT_3AT
hdf5file_3AT <- rhdf5::H5Fopen(hd_file_3AT) # filehandle for the h5 file

## WT_CHX
WT_CHX <- "G-Sc_2014/output/WTCHX/WTCHX.h5"
hd_file_CHX <- WT_CHX
hdf5file_CHX <- rhdf5::H5Fopen(hd_file_CHX) # filehandle for the h5 file

asite_disp_length_file <- "G-Sc_2014/input/asite_disp_length_yeast_standard.txt"

# asite_disp_length <- readr::read_tsv(asite_disp_length_file,
#                                      comment = "#"
# )


#Initial set ups
dataset_G2014 <- "G-Sc_2014"
test_orfs <- c("YCR012W","YEL009C","YOR303W","YOL130W","YGR094W", "YML065W")
genes <- c("YPR036W-A", "YEL009C", "YOR303W", "YGR094W", "YOL130W", "YOR061W", "YAL040C", "YKL109W", "YIL144W", "YGR037C", "YJL106W", "YBR160W", "YDR172W", "YML065W", "YBR257W", "YOL104C", "YJR094C","YCR012W")
gene_names <- rhdf5::h5ls(hdf5file_none, recursive = 1)$name

nnt_gene<- 50
startpos <-250
startlen <- 10
temporary_length <- 20 #for zooming in bit
min_read_count <- 10
orf_gff_file <- "G-Sc_2014/input/yeast_CDS_w_250utrs.gff3"
# Edward_gff <- "yeastutrgff/out/transcriptcentric_abundant_full-ORF_ypd_plus_other_fixed_UTR_length_transcripts.gff"

#Function to create a gff table  
readGFFAsDf <- purrr::compose(
  rtracklayer::readGFFAsGRanges,
  data.frame, 
  as_tibble,
  .dir = "forward" # functions called from left to right
)

gff_df <- readGFFAsDf(orf_gff_file) 
# gff_Edward <- readGFFAsDf(Edward_gff) 

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
Get5UTRstart <- function(gene, gffdf) {
  gffdf %>% 
    dplyr::filter(Name == gene, type=="UTR5" ) %>% 
    dplyr::pull(start)
}

# Get5UTRstart <- function(gene, gffdf) {
#   gffdf %>%
#     dplyr::filter(str_detect(seqnames, gene), type=="five_prime_UTR" )  %>%
#     dplyr::pull(start)
# }

# Takes value for 3' end of the table (includes whole 5'UTR +50 nt of CDS)
CDS3_end <- function(gene, gffdf) {
  gffdf %>% 
    dplyr::filter(Name == gene, type=="UTR5") %>% 
    dplyr::pull(end) + nnt_gene
} 

# CDS3_end <- function(gene, gffdf) {
#   gffdf %>%
#     dplyr::filter(str_detect(seqnames, gene), type == "five_prime_UTR" ) %>%
#     dplyr::pull(end) + nnt_gene
# }

# Takes value for length of 5'UTR  
UTR5_length <- function(gene, gffdf) {
  gffdf %>% 
    dplyr::filter(Name == gene, type=="UTR5" ) %>% 
    dplyr::pull(width)
} 


# UTR5_length <- function(gene, gffdf) {
#   gffdf %>%
#     dplyr::filter(str_detect(seqnames, gene), type == "five_prime_UTR" ) %>%
#     dplyr::pull(width)
# }


# Creates a matrix with number of columns that start at 5' UTR and finish at 50'th nt of the CDS
GetGeneDatamatrix5UTR <- function(gene, dataset, hdf5file, gffdf, nnt_gene) {
  data_mat_all <- GetGeneDatamatrix(gene, dataset, hdf5file)
  n_left5 <- Get5UTRstart(gene, gffdf) # column to start from (5'end)
  n_right3 <-CDS3_end(gene, gffdf)  # column to end with (3'end)
  data_mat_5start <- data_mat_all[, n_left5 : n_right3]
  # data_mat_5start <- tibble::as_tibble(data_mat_5start, .name_repair="minimal")
  return(data_mat_5start)
}
# 
# # GetGeneDatamatrix5UTR(test_orfs[1], dataset = dataset, hdf5file_none, gff_df, nnt_gene = nnt_gene)
# 
# GetPosnCountOutput <- function(x){
#   output_thing <- x %>%
#     gather(-read_length, key="Position", value = "Counts") %>%
#     group_by(Position) %>%
#     summarise(Counts=sum(Counts)) %>%
#     arrange(as.integer(Position))
#   return(output_thing)
# }
# 
# TidyDatamatrix <- function(x, startpos = 1, startlen = 1, gene) {
#   # CHECK startpos/off-by-one
#   positions <- startpos:(startpos + ncol(x) - 1)
#   readlengths <- startlen:(startlen + nrow(x) - 1)
#   x %>%
#     magrittr::set_colnames(positions) %>%
#     tibble::as_tibble() %>%
#     dplyr::mutate(ReadLen = readlengths) %>%
#     tidyr::gather(-ReadLen, key = "Pos", value = "Counts", convert = FALSE) %>%
#     dplyr::mutate(Pos = as.integer(Pos), Counts = as.integer(Counts)) %>%
#     dplyr::group_by(Pos) %>%
#     dplyr::summarise(Counts=sum(Counts))
# }
# 
# # final function from raw data processing to data visualization
# GetGeneDatamatrix5UTR_non_A <- function(gene, dataset, hdf5file, gffdf, nnt_gene) {
#   data_mat_all <- GetGeneDatamatrix(gene, dataset, hdf5file)
#   n_left5 <- Get5UTRstart(gene, gffdf) # column to start from (5'end)
#   n_right3 <- UTR5_length(gene, gffdf) + nnt_gene # column to end with (3'end)
#   data_mat_5start <- data_mat_all[, n_left5 : n_right3]
#   data_mat_5start <- tibble::as_tibble(data_mat_5start, .name_repair="minimal")
#   return(data_mat_5start)
# }
# 
# UTR5_table <- function(gene) {
#  lapply(gene,
#          function(gene)
#            GetGeneDatamatrix5UTR_non_A(gene,
#                                  dataset,
#                                  hdf5file,
#                                   gffdf =gff_df,
#                                  nnt_gene = 50)
#   ) %>%
#     Reduce("+", .) %>% # sums the list of data matrices
#     TidyDatamatrix(startpos = -250, startlen = 10) %>%
#  return()
# 
# }
# # 
# abc2 <- UTR5_table(test_orfs[2])
# 
# # the final function i'd use :) iteration of UTR5_table
# final_function_table <- function(genes){
#  table <-purrr::map(genes, UTR5_table)
#  names(table) <- genes
#  return(table)
# 
# }
# output_orfs<-final_function_table(test_orfs) #it works just fine
# # 
# # ##meta-analysis:
# # 
# # meta_5genes_none <- UTR5_table(genes)
# #   ##just do the plotting here
# # 
# # meta_5genes_CHX <- UTR5_table(genes)
# # 
# plotting_meta_analysis<- function(input_data) {
#   # text_AUG <- textGrob("AUG", gp=gpar(fontsize=13, fontface="bold")) #to place text annotation
# 
#   #the actual plotting
#   plotted_UTR <- ggplot(input_data) +
#     geom_density(aes(x=Pos, y=Counts), stat="identity") +
#     scale_x_continuous(limits = c(-250,50), expand = c(0, 0)) +
#     scale_y_continuous(expand = c(0,0)) +
#     labs(y= "Read count", x = "Position") +
#      ggtitle(paste0("Ribosome footprint density of all genes: WT 3-AT")) +
#     # coord_cartesian(clip = "off") +
#     # annotation_custom(text_AUG,xmin=0,xmax=0,ymin=0,ymax=5) +
#     theme_classic() %>%
#     return()
# }
# #for WT_none
# meta_5genes_plot_none <- plotting_meta_analysis(meta_5genes_none)
# ggsave("5UTR_notAmapped", device = "jpg")
# 
# #for WT_CHX
# meta_5genes_plot_CHX <- plotting_meta_analysis(meta_5genes_CHX)
# ggsave("5UTR_notAmapped_CHX", device = "jpg")
# ##I should still be using these functions as they will allow me to compare the A- and non-A site mappng!

#########################A site displacement slot#########################

###
## scaling factor:
GetGeneReadsTotal <- function(gene, dataset, hdf5file) {
  rhdf5::H5Aread(rhdf5::H5Aopen(rhdf5::H5Gopen(hdf5file, paste0("/", gene, "/", dataset, "/reads")), "reads_total"))
}

scaling_factor_none <<- sapply(gene_names, GetGeneReadsTotal, dataset_G2014, hdf5file_CHX) %>% sum()/ 10^6

scaling_factor_CHX <<- sapply(gene_names, GetGeneReadsTotal, dataset_G2014, hdf5file_none) %>% sum()/ 10^6

scaling_factor_3AT <<- sapply(gene_names, GetGeneReadsTotal, dataset_G2014, hdf5file_3AT) %>% sum()/ 10^6

####

CalcAsiteFixedOneLength <- function(reads_pos_length, min_read_length,
                                    read_length, asite_disp) {
  # Calculate read A-site using a fixed displacement for a single read length
  length_row_choose <- read_length - min_read_length + 1
  reads_pos_length[length_row_choose, ] %>%
    dplyr::lag(n = asite_disp, default = 0)
}

# CalcAsiteFixedOneLength(reads_pos_length = GetGeneDatamatrix(gene = test_orfs[1],dataset = dataset,hdf5file = hdf5file), min_read_length = 10, read_length = 28, asite_disp = 15)
# > str(something)
# num [1:1751] 0 0 0 0 0 0 0 0 0 0 ...

###
CalcAsiteFixed <- function(reads_pos_length, min_read_length,
                           asite_disp_length,
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

# CalcAsiteFixed(reads_pos_length = GetGeneDatamatrix(gene = test_orfs[1], dataset = dataset, hdf5file = hdf5file), min_read_length, min_read_length = 10, colsum_out = TRUE)

###
SnapToCodon <- function(x, left, right, snapdisp=0L) {
  # snap nucleotide-aligned reads to codon position
  #   x:     vector
  #   left:  integer for starting position, frame 0
  #   right: integer for ending position
  #   snapdisp: integer any additional displacement in the snapping
  RcppRoll::roll_suml(x[(left:right) + snapdisp], n=1L, by=1L, fill = NULL)
  
}  

#
###
# GetGeneCodonPosReads1dsnap <- function(gene, dataset, hdf5file, left, right,
#                                        min_read_length,
#                                        asite_disp_length = data.frame(
#                                          read_length = c(28, 29, 30),
#                                          asite_disp = c(15, 15, 15)
#                                        ),
#                                        snapdisp=0L) {
#   left = Get5UTRstart(gene, gff_df)
#   right = CDS3_end(gene, gff_df)
# 
#   reads_pos_length <- GetGeneDatamatrix(gene, dataset, hdf5file)
#   reads_asitepos <- CalcAsiteFixed(
#     reads_pos_length, min_read_length,
#     asite_disp_length
#   )
# 
#   Counts <- SnapToCodon(reads_asitepos,left,right,snapdisp)
#   Pos <- UTR5_table(gene) %>%
#     select( Pos)
# 
#   cbind(Pos, Counts) %>% as_tibble()
# 
# }

### oh noooo, that was suppsed to be changed!
# GetPosition <- function(x, startpos, endposition)  {
#   posStart <- (startpos - endposition)
#   
#   positions <- posStart:(posStart + ncol(x) - 1) 
#   positions %>% 
#     tibble::enframe(name = NULL) %>%
#     set_colnames("Pos")
# }  

GetPosition <- function(x, startpos = 1)  {
  positions <- startpos:(startpos + ncol(x) - 1) 
  positions %>% 
    tibble::as_tibble() %>%
    magrittr::set_colnames("Pos")
}  

###
GetGeneCodonPosReads1dsnap <- function(gene, dataset, hdf5file, gffdf,
                                       nnt_gene, min_read_length, asite_disp_length,  snapdisp = 0L, scaling_factor) {
  
  left  = Get5UTRstart(gene, gffdf)
  right = CDS3_end(gene, gffdf) 
  # end = UTR5_length(gene, gffdf)
  
  # will have to call x and gff_df the same!
  reads_pos_length <- GetGeneDatamatrix5UTR(gene, dataset,
                                            hdf5file, gffdf,nnt_gene)
  
  
  reads_asitepos <- CalcAsiteFixed(reads_pos_length,
                                   min_read_length,
                                   asite_disp_length)  
  
  Counts <- SnapToCodon(reads_asitepos,
                        left = left,
                        right = right,
                        snapdisp) 
  
  Pos <- reads_pos_length %>%
    GetPosition(startpos = -250)
  
  final_tibble <- base::cbind(Pos, Counts) %>% 
    tibble::as_tibble() 
  
  norm_tibble <<- final_tibble %>%
    mutate(., Normalized_Counts = apply(.[,"Counts"], 1, function(x) sum(x)/scaling_factor )) %>%
    select(Pos, Normalized_Counts) %>%
    set_colnames(c( "Pos", "Counts"))
  return(norm_tibble)
  
}


output1 <-GetGeneCodonPosReads1dsnap(gene = test_orfs[1], dataset = dataset_G2014, hdf5file = hdf5file_CHX, gffdf = gff_df, nnt_gene, min_read_length = 10, asite_disp_length = asite_disp_length, snapdisp = 0L, scaling_factor = scaling_factor_CHX) 

output2 <- GetGeneCodonPosReads1dsnap(gene = test_orfs[2], dataset = dataset_G2014, hdf5file = hdf5file_CHX, gffdf = gff_df, nnt_gene, min_read_length = 10, asite_disp_length = asite_disp_length, snapdisp = 0L, scaling_factor = scaling_factor_CHX) 

combined<-full_join(output1, output2, by = "Pos")
combined %>%
  mutate(., Counts = Counts.x + Counts.y)


#  map(test_orfs, ~GetGeneCodonPosReads1dsnap(
#   .,
#   dataset,
#   hdf5file,
#   gff_df,
#   nnt_gene = 50,
#   min_read_length = 10,
#   asite_disp_length = asite_disp_length ))
# names(Counts_Asite_mapped) <- test_orfs

#for multiple genes, the output are singular tibbles
A_mapped_genes <- function(gene,
                           dataset,
                           hdf5file,
                           gffdf,
                           min_read_length, 
                           asite_disp_length,
                           scaling_factor) {
  
  
  output<- purrr::map(gene,
                      GetGeneCodonPosReads1dsnap,
                      dataset,
                      hdf5file,
                      gffdf,
                      nnt_gene,
                      min_read_length,
                      asite_disp_length,
                      scaling_factor = scaling_factor_3AT)
  
  names(output) <- gene
  
  # filter(.data = output, Pos >= -250, Pos < 0) 
  # group_by(Pos)%>%
  #   summarise(sum_read_counts_uAUG=sum(Counts))
  # 
  return(output)
}
######REALLY IMPORTAT...

########### trying to combine multiple datasets, maybe i should do that at the matrix level?
# mapped_A_genes_CHX <- A_mapped_genes(test_orfs, dataset = dataset_G2014, hdf5file = hdf5file_CHX, gffdf = gff_df, min_read_length = 10, asite_disp_length = asite_disp_length)
# 
# mapped_A_genes_CHX <-bind_rows(mapped_A_genes_CHX,.id = "Gene")
# 
# mapped_A_genes_none <- A_mapped_genes(test_orfs, dataset = dataset_G2014, hdf5file = hdf5file_none, gffdf = gff_df, min_read_length = 10, asite_disp_length = asite_disp_length)
# 
# mapped_A_genes_none <-bind_rows(mapped_A_genes_none,.id = "Gene")
# 
# joining_datasets <-function(dataset1, dataset2, gene){
#   data1 <- filter(dataset1, Gene == gene)
#   data2 <- filter(dataset2, Gene == gene)
#   
#   tibble_created <-full_join(data1, data2, by = "Pos") %>%
#     mutate(., Counts = Counts.x + Counts.y) %>%
#     select( Pos, Counts) 
#   names(tibble_created) <- gene
#   return(tibble_created)
# }
# 
# joining_datasets(mapped_A_genes_CHX, mapped_A_genes_none, "YGR094W")
# 
# map(test_orfs, joining_datasets, mapped_A_genes_CHX, mapped_A_genes_none)
# 
# joining_datasets <-function(dataset1, dataset2, gene){
#    dataset1 <- dataset1[[gene]]
#    dataset2 <- dataset2[[gene]]
#   
#   full_join(dataset1, dataset2, by = "Pos") %>%
#     return()
# }

###############


# normalization_AUG_UTR <- function(mapped_dataset){
#   
#   #normalize read counts mapped to each gene for it's length (5'UTR + 50nt of CDS)
#   norm_plan1 <- mapped_dataset %>%
#     bind_rows( .id = "Gene") %>%
#     filter( Pos >= -250, Pos < -5) %>%
#     as_tibble() %>%
#     dplyr::group_by(Gene) %>%
#     dplyr::summarise(All_read_counts = sum(Counts)/260) 
#   
#   #find the scaling factor;sum(read counts from all genes)/10^6
#   scaling_factor1 <-norm_plan1 %>%
#     summarise(sum(All_read_counts)/10^6) %>%
#     pull
#   
#   #divide the normalized read counts by the scaling factor (normalize for seq depth)
#   Normalized_Reads_5UTR <-norm_plan1 %>%
#     group_by(Gene) %>%
#     summarise(Normalized_Reads_5UTR = (All_read_counts/scaling_factor1)/245) %>%
#     return()
#   ## divide by 245 to get the density of reads in the 5'UTR
#   
#   norm_plan2 <-mapped_dataset %>%
#     bind_rows( .id = "Gene") %>%
#     filter( Pos >= -5, Pos <= 10) %>%
#     as_tibble() %>%
#     dplyr::group_by(Gene) %>%
#     dplyr::summarise(All_read_counts = sum(Counts)/260) 
#   
#   #find the scaling factor;sum(read counts from all genes)/10^6
#   scaling_factor2 <-norm_plan2 %>%
#     summarise(sum(All_read_counts)/10^6) %>%
#     pull
#   
#   #divide the normalized read counts by the scaling factor (normalize for seq depth)
#   Normalized_Reads_AUG <-norm_plan2 %>%
#     group_by(Gene) %>%
#     summarise(Normalized_Reads_AUG = (All_read_counts/scaling_factor2/15)) %>%
#     return()  
#   ## divide by 15 to get the density of reads in the AUG
#   
#   full_join(Normalized_Reads_AUG, Normalized_Reads_5UTR, by = "Gene")
#   
# }
#############################################################################
#CanVsNon - canonical vs non-canonical
CanVsNon <- function(dataset) {
  
  UTR <- bind_rows(dataset, .id = "Gene") %>%
    dplyr::filter( Pos >= -250, Pos < -5) %>%
    dplyr::group_by(Gene) %>%
    dplyr::summarise(Counts_UTR5 = sum(Counts)/(length(-250:-5)-1))
  
  AUG <- bind_rows(dataset, .id = "Gene") %>%
    dplyr::filter( Pos >= -5, Pos < 10) %>%
    dplyr::group_by(Gene) %>%
    dplyr::summarise(Counts_AUG = sum(Counts)/(length(-5:10)-1))
  
  UTR_AUG <-dplyr::full_join(UTR, AUG, by = "Gene") 
  
}

all_genes <- A_mapped_genes(test_orfs, dataset = dataset_G2014, hdf5file = hdf5file_3AT, gffdf = gff_df, min_read_length = 10,asite_disp_length)

all_genes_plot <- CanVsNon(all_genes) 

###a scatter plot showing the UTR/AUG ratio
ggplot(all_genes_plot, aes(x = Counts_AUG, y = Counts_UTR5)) +
  geom_point() +
  scale_y_continuous(expand = c(0,0)) +
  scale_y_log10() +
  geom_text_repel(aes(label = Gene))


#####
sevR1 <-A_mapped_genes(gene = test_orfs[1:2], dataset = dataset_G2014, hdf5file = hdf5file_none, gffdf = gff_df, min_read_length = 10, asite_disp_length)

sevR1 <- CanVsNon(sevR1)

sevR2 <-A_mapped_genes(gene = test_orfs[1:2], dataset = dataset_G2014, hdf5file = hdf5file_CHX, gffdf = gff_df, min_read_length = 10,asite_disp_length = asite_disp_length)

sevR2 <- CanVsNon(sevR2)

sevR3 <-A_mapped_genes(gene = test_orfs[1:2], dataset = dataset_G2014, hdf5file = hdf5file_3AT, gffdf = gff_df, min_read_length = 10,asite_disp_length)

sevR3 <- CanVsNon(sevR3)

comparing_UTR_test <- full_join(sevR1, sevR2, by = "Gene") %>%
  full_join(sevR3, by = "Gene") %>%
  set_colnames(c("Gene", "UTR_PRE_ENTRY_1", "AUG_PRE_ENTRY_1", "UTR_META_I_1","AUG_META_I_1", "UTR_META_II_1", "AUG_META_II_1" )) 

comparing_UTR_test_1 <-comparing_UTR_test%>%
  select("Gene","UTR_PRE_ENTRY_1","UTR_META_I_1","UTR_META_II_1") %>%
  gather(key = "Condition", value = "Reads", -Gene)

ggplot(comparing_UTR_test_1, aes(x = Condition,group = 1, y= Reads, color = Gene)) +
  geom_point() +
  geom_line() +
  geom_text_repel(aes(label = round(Reads, digits = 3)))


# ####
# ### calculate read counts at each position for all genes 
All_genesAmapped <- function(gene, dataset, hdf5file, gffdf,
                             min_read_length, asite_disp_length,
                             scaling_factor) {
  
  output<- purrr::map(gene,
                      GetGeneCodonPosReads1dsnap,
                      dataset,
                      hdf5file,
                      gffdf,
                      nnt_gene,
                      min_read_length,
                      asite_disp_length,
                      scaling_factor = scaling_factor_3AT)
  
  names(output) <- gene
  
  Counts <- output %>%
    Reduce("+", .) %>%
    dplyr::select(Counts)
  
  Pos <- output$YEL009C$Pos ### !!!
  
  Counts_Asite_mapped_all <-base::cbind(Pos, Counts) %>%
    tibble::as_tibble(.name_repair = "minimal")
  
  return(Counts_Asite_mapped_all)
  
}

# 
try <- All_genesAmapped(gene = test_orfs, dataset = dataset_G2014, hdf5file = hdf5file_CHX, gffdf = gff_df, min_read_length = 10, asite_disp_length = asite_disp_length)

try2 <- All_genesAmapped(gene = test_orfs, dataset = dataset_G2014, hdf5file = hdf5file_none, gffdf = gff_df, min_read_length = 10, asite_disp_length = asite_disp_length)

join_reduced_data <-function(dataset1, dataset2){
  full_join(dataset1,dataset2, by = "Pos") %>%
    mutate(., Counts = Counts.x + Counts.y) %>%
    select(Pos, Counts)%>%
    return()
} 
## this one works which is really good!


# 
# try <- normalization_each_gene(mapped_dataset = mapped_A_genes)
# 
# 
# ############################################################################################                                               ##plotting##
# 
plotting_5UTR<- function(input_data) {
  # text_AUG <- textGrob("AUG", gp=gpar(fontsize=13, fontface="bold")) #to place text annotation
  
  #the actual plotting
  plotted_UTR <- ggplot2::ggplot(input_data) +
    geom_density(aes(x=Pos, y=Counts), stat="identity") +
    scale_x_continuous(limits = c(-250, 50), expand = c(0, 0)) +
    scale_y_continuous(expand = c(0,0)) +
    labs( x = "Position", y = "Footprint density (CPM)") +
    # ggtitle(paste0("Ribosome footprint density of ", input_data)) +
    # coord_cartesian(clip = "off") +
    # annotation_custom(text_AUG,xmin=0,xmax=0,ymin=0,ymax=5) +
    theme_classic()
}

plotting_5UTR_with_title<- function(input_data) {
  # text_AUG <- textGrob("AUG", gp=gpar(fontsize=13, fontface="bold")) #to place text annotation
  
  #the actual plotting
  plotted_UTR <- ggplot2::ggplot(input_data) +
    geom_density(aes(x=Pos, y=Counts), stat="identity") +
    scale_x_continuous(limits = c(-250,50), expand = c(0, 0)) +
    scale_y_continuous(expand = c(0,0)) +
    labs( x = "Position") +
    ggtitle(paste0("Ribosome footprint density of WT_3AT")) +
    # coord_cartesian(clip = "off") +
    # annotation_custom(text_AUG,xmin=0,xmax=0,ymin=0,ymax=5) +
    theme_classic()
}

##multiple plots##

plotting_multiple <- function(input_data, gene_names) {
  purrr :: map(input_data, plotting_5UTR) %>%
    ggarrange(plotlist = ., labels = gene_names) %>%
    return()
  #ggsave(file.path(paste0("Plot of ", gene)), device = "jpg")
}
# 
# ####################################################################################
# ##2 in one##
# 
# #this is what i'd like to achieve for the final A-site mapping function 
# # choosing_sample <- function(dataset, genes) {
# # hd_file <- dataset
# # hdf5file <- rhdf5::H5Fopen(hd_file) # filehandle for the h5 file
# # 
# # final_function_table(genes) %>%
# # return()
# # }
# 
# # choosing_sample(WT_none, genes = test_orfs)
# 
# ##WT none
# # WT_none <- "G-Sc_2014/output/WTnone/WTnone.h5"
# # hd_file <- WT_none
# # hdf5file <- rhdf5::H5Fopen(hd_file) # filehandle for the h5 file
# # output_none <-final_function_table(test_orfs)
# # 
# # ## WT_3AT
# # WT_3AT <- "G-Sc_2014/output/WT3AT/WT3AT.h5"
# # hd_file <- WT_3AT
# # hdf5file <- rhdf5::H5Fopen(hd_file) # filehandle for the h5 file
# # output_3AT <-final_function_table(test_orfs)
# # 
# # ## WT_CHX
# # WT_CHX <- "G-Sc_2014/output/WTCHX/WTCHX.h5"
# # hd_file <- WT_CHX
# # hdf5file <- rhdf5::H5Fopen(hd_file) # filehandle for the h5 file
# # output_CHX <-final_function_table(test_orfs)
# ######
# 
two_in_one_plots <- function(sample1, sample2, gene, cond1_name, cond2_name) {
  ggplot2::ggplot(sample1,aes(x=Pos, y=Counts)) +
    geom_density(data = sample1, stat="identity", aes(color = "VEG")) +
    geom_density(data = sample2, stat="identity", aes(color = "MEIOSIS")) +
    scale_x_continuous(limits = c(-250,50), expand = c(0, 0)) +
    scale_y_continuous(expand = c(0,0)) +
    labs(y= "Read count", x = "Position") +
    theme_classic() +
    ggtitle(paste0(gene, ":", cond1_name, " vs ", cond2_name)) +
    scale_color_manual(values = c(VEG = "red", MEIOSIS = "blue"))
}

three_in_one_plots <- function(sample1, sample2, sample3, gene) {
  ggplot2::ggplot(sample1,aes(x=Pos, y=Counts)) +
    geom_density(data = sample1, stat="identity", aes(color = "VEG")) +
    geom_density(data = sample2, stat="identity", aes(color = "ANAPH_I")) +
    geom_density(data = sample3, stat="identity", aes(color = "SPORE")) +
    scale_x_continuous(limits = c(-250,50), expand = c(0, 0)) +
    scale_y_continuous(expand = c(0,0)) +
    labs(y= "Read count", x = "Position") +
    theme_classic() +
    scale_color_manual(values = c(VEG = "red", ANAPH_I = "blue", SPORE = "green" ))
}


# #### that's the kind of stuff i'd like to achieve 
# 
# # WT_CHX_3AT <- two_in_one_plots(output_3AT[[1]], output_CHX[[1]], "YCR012W", "WT 3AT", "WT CHX") 
# # they're still lists and we want tibbles, otherwise it doesnt work
# ##I have no idea how to iterate this one
# 
# ####################################################################################
# ### Creating a table for multiple positions for multiple genes
# 
#creating adjusted regions to sum read counts
sum_uAUG <- function(datatibble, uAUG_start, uAUG_end){
  datatibble %>%
    dplyr::filter(Pos >= uAUG_start, Pos <= uAUG_end) %>%
    dplyr::summarise(sum_read_counts_uAUG=sum(Counts)/(length(uAUG_start:uAUG_end)-1)) 
  # dplyr::pull
}

#for multiple genes
sum_uAUG_multiple<- function(genes, uAUG_start, uAUG_end) {
  #function 
  sapply(genes, sum_uAUG, uAUG_start, uAUG_end) %>%
    return()
  #sapply jest lepsze bo nie tworzy listy ale jakieś problemy są 
  # set_rownames(., paste0(gene_names)) %>%
  #names of columns and rows
}
# 
# #####CDS VS 5'UTR:
# 
CDS_5UTR <- function(genes, gene_names) {
  table_genes <-as.character((paste0(gene_names)))  #why does it make it a factor?
  
  
  region1 <-sum_uAUG_multiple(genes, uAUG_start = -250, uAUG_end = -5)
  region2 <- sum_uAUG_multiple(genes, uAUG_start = -5, uAUG_end = 10)
  
  
  base::cbind(table_genes, region1, region2) %>%
    as_tibble() %>%
    set_colnames(c("genes","5'UTR", "AUG")) %>%
    # tidyr::gather( key = "region", value = "count", -genes) %>%
    return()
}

#required for efficiency barplot
together <-CDS_5UTR(genes = META_I_1_sev,gene_names = genes) %>%
  tidyr::gather( key = "region", value = "count", -genes)
together$count <- as.numeric(together$count)
# 
# 
# #########
efficiency_barplot <- function(data) {
  positions <- c( "5'UTR", "AUG")
  
  ggplot(data,aes(x=region, y = count, fill = region )) +
    geom_col(show.legend = FALSE) +
    theme_classic() +
    labs(y= "Footprint density", x = "Region (mRNA)") +
    # ggtitle("quantification across all genes") +
    scale_y_continuous(expand = c(0,0)) +
    scale_x_discrete(limits = positions) 
  
  ## that gives an overall view of all read counts at each region
}

efficiency_barplot(WT_NONE_sev_UTR_AUG_plot)


# 
efficiency_barplot_each_gene <- function(data, gene){
  positions <- c( "5'UTR", "AUG")
  value <-filter(data, genes == gene)
  
  ggplot(value, aes(x = region, y = count, fill = region)) +
    geom_col(show.legend = FALSE) +
    theme_classic() +
    labs(y= "Footprint density (CPMB)", x = "Region (mRNA)") +
    ggtitle(paste0(gene)) +
    scale_y_continuous(expand = c(0,0)) +
    scale_x_discrete(limits = positions)
}

efficiency_barplot_each_gene(together, "YJR094C")

# map(together, efficiency_barplot_each_gene, gene = genes )
# # Error in UseMethod("filter_") : 
# #   no applicable method for 'filter_' applied to an object of class "list"
# 
# #########
# #creating table for all genes together 
# # together_all_genes <-CDS_5UTR(genes = mapped_A_genes,gene_names = test_orfs)
# # together_all_genes$`5'UTR` <- as.numeric(together_all_genes$`5'UTR`)
# # together_all_genes$AUG <- as.numeric(together_all_genes$AUG)
# # colSums(together_all_genes[,-1]) 
# # 5'UTR   AUG 
# #   177    36
# 
# ##################
# #0:-15 nt
# #0:-60 nt
# #60:-120 nt
# #0:-250 nt
all_regions_together <- function(genes, gene_names) {
  table_genes <-as.character((paste0(gene_names)))  #why does it make it a factor?
  
  region1 <-sum_uAUG_multiple(genes, uAUG_start = -250, uAUG_end = 0)
  region2 <- sum_uAUG_multiple(genes, uAUG_start = -5, uAUG_end = 10)
  region3 <-sum_uAUG_multiple(genes, uAUG_start = -15, uAUG_end = 0)
  region4 <-sum_uAUG_multiple(genes, uAUG_start = -60, uAUG_end = 0)
  region5 <-sum_uAUG_multiple(genes, uAUG_start = -120, uAUG_end = -60)
  region6  <-sum_uAUG_multiple(genes, uAUG_start = -250, uAUG_end = -120)
  
  # region5 <- sum_uAUG_multiple(genes, gene_names, uAUG_start = -15, uAUG_end = 10) %>% unlist(as.numeric())
  
  base::cbind(table_genes, region1, region2, region3, region4, region5) %>%
    as_tibble() %>%
    magrittr::set_colnames(c("genes","5'UTR", "AUG", "-15:-0", "-60:0", "-120:-60")) %>%
    tidyr::gather( key = "region", value = "read", -genes) %>%
    return()
}
# 
# #example: 
# 
together <- A_mapped_genes(test_orfs, dataset = dataset_G2014, hdf5file = hdf5file_CHX, gffdf = gff_df, min_read_length = 10, asite_disp_length = asite_disp_length) %>%
  all_regions_together(gene_names = test_orfs)


WT_NONE_sev_regions$read <- as.numeric(WT_NONE_sev_regions$read)

plotted_barplots(WT_NONE_sev_regions)

# # ## that gives an overall view of all read counts at each region 
plotted_barplots <- function(data) {
  positions <- c( "5'UTR", "AUG", "-15:-0", "-60:0", "-120:-60", "-250:-60")
  
  ggplot(data,aes(x=region, y = read, fill = region)) +
    geom_col(show.legend = FALSE) +
    theme_classic() +
    labs(y= "Footprint density (CPMB)", x = "Region (mRNA)") +
    ggtitle("quantification across all genes") +
    scale_y_continuous(expand = c(0,0)) +
    scale_x_discrete(limits = positions)
  
  ## that gives an overall view of all read counts at each region
}
# 
# ## Individual for each gene
plotted_each_barplot <- function(data, gene){
  positions <- c( "5'UTR", "AUG", "-15:-0", "-60:0", "-120:-60","-250:-60")
  value <-filter(data, genes == gene)
  
  ggplot(value,aes(x = region, y = read, fill = region)) +
    geom_col(show.legend = FALSE) +
    theme_classic() +
    labs(y= "Footprint density (CPMB)", x = "Region (mRNA)") +
    ggtitle(paste0(gene)) +
    scale_y_continuous(expand = c(0,0)) +
    scale_x_discrete(limits = positions)
}



# ##########################################################################################
# 
### Comparison of 5'UTR for all genes between 3 conditions
# requires PRE_ENTRY_2_sev_5UTR as the output

UTR5_3_conditions_all <-function(condition1_UTR, condition2_UTR, condition3_UTR, name1, name2, name3){
  condition1 <- sum(condition1_UTR$Counts_UTR5)
  condition2 <- sum(condition2_UTR$Counts_UTR5)
  condition3 <- sum(condition3_UTR$Counts_UTR5)
  
  Condition <- paste0(c(name1, name2, name3))
  # Condition <- c("WT NONE", "WT CHX", "WT 3AT")
  for_plot_UTR <- rbind(condition1, condition2, condition3) %>%
    as_tibble %>%
    set_colnames("UTR5") %>%
    cbind(Condition)
  
  plot_UTR <-ggplot(for_plot_UTR) +
    geom_bar(aes(x=Condition, y= UTR5, fill = Condition), stat="identity", show.legend = FALSE ) +
    # scale_x_discrete(labels=c("WT_3AT","WT_CHX", "WT_NONE")) +
    scale_y_continuous(expand = c(0,0)) +
    theme_classic() +
    labs(x= "Stage", y = "5'UTR Footprint density (CPMB)")
  return(plot_UTR)
}

##### for all genes
all_genes_UTR_AUG <-function(condition1_UTR){
  UTR <- sum(condition1_UTR$Counts_UTR5)
  AUG <- sum(condition1_UTR$Counts_AUG)
  
  Condition <- c("5 UTR", "AUG")
  for_plot_UTR <- rbind(UTR, AUG) %>%
    as_tibble %>%
    set_colnames("UTR5") %>%
    cbind(Condition)
  
  plot_UTR <-ggplot(for_plot_UTR) +
    geom_bar(aes(x=Condition, y= UTR5, fill = Condition), stat="identity", show.legend = FALSE ) +
    # scale_x_discrete(labels=c("WT_3AT","WT_CHX", "WT_NONE")) +
    scale_y_continuous(expand = c(0,0)) +
    theme_classic() +
    labs(x= "Region", y = "5'UTR Footprint density (CPMB)")
  return(plot_UTR)
}
all_genes_UTR_AUG(VEG_ALL_GENES_5UTR)


### Comparison of 5'UTR for 1 gene between 3 conditions
# requires PRE_ENTRY_2_sev_5UTR as the output, also gene must be entre parentheses

compared_UTR_single<-function(condition1, condition2, condition3, gene){
  
  search_for_5UTR <- function(condition, gene){
    condition %>%
      filter(Gene == gene ) %>%
      select(Counts_UTR5) %>%
      pull
  }
  
  UTR1 <-search_for_5UTR(condition1, gene)
  
  UTR2 <-search_for_5UTR(condition2, gene)
  
  UTR3 <-search_for_5UTR(condition3, gene)
  
  Condition <- c("WT NONE", "WT CHX", "WT 3AT")
  
  for_plot <- rbind(UTR1,UTR2,UTR3) %>%
    set_colnames("UTR5") %>%
    cbind(Condition) %>%
    as_tibble()
  for_plot$UTR5 <-as.numeric(for_plot$UTR5)
  
  plot_UTR <-ggplot(for_plot) +
    geom_bar(aes(x=Condition, y= UTR5, fill = Condition), stat="identity",show.legend = FALSE ) +
    # scale_x_discrete(labels=c("META_I_1", "META_II_1", "PRE_ENTRY_1")) +
    scale_y_continuous(expand = c(0,0)) +
    theme_classic() +
    labs(x= "Condition", y = "5'UTR Footprint density (CPMB)") +
    ggtitle(paste0(gene))
  
  return(plot_UTR)
  
}
# 
# 
# ##########################################################################################
#                                      ######FASTQ######
# fq.file <- "/Users/Ania/Desktop/Szkoła/4th year/Dissertation/gene_graphs/Ania_Riboviz/G-Sc_2014/input/yeast_CDS_w_250utrs.fa"
# 
# fastq_sequence <- Fastqfile(fq.file)
#   #not entirely what i want, a bit confusing: 
# # A tibble: 37,762 x 3
# # Header                        Sequence                      Quality                   
# # <chr>                         <chr>                         <chr>                       
# #   1 >YAL068C                      ACCTATGAAAGATTTATGATTCGTTCAG… GCTGCTTCAACTATATGCCTTTGAGAAT…
# # 2 AACAAATACAATGGTCAAATTAACTTCA… TAGCTCAATCTGACGAAAGAGTCAACTT… CACCATGTTGACCGGTATTGCTCCAGAC…
# # 3 CAGCCATCTCCAGTGCTCTATCCAAGGA… TTCCATAGAAATTGAAAATTAACGAACA… AAAGAAACTTCTACACTATTGTAGAAAA…
# # 4 >YAL067W-A                    ATATTCTCAAAGGCATATAGTTGAAGCA… TTCTGAACGAAT
# 
# #################### FULL BRAR SECTION ########################################
# #### Brar et al. 2012 data
# choosing_sample <- function(hdfile, genes, dataset, gffdf, min_read_length) {
#   
# hd_file <- hdfile
# hdf5file <- rhdf5::H5Fopen(hd_file) # filehandle for the h5 file
# 
# All_genesAmapped(genes, dataset, hdf5file, gffdf, min_read_length) %>%
# return()
# }
# 
# choosing_sample(hd_file_PRE1, test_orfs, dataset = dataset_Brar, gff_df, min_read_length = 10) #works!! still gotta improve the file part
# 
# 
# dataset_Brar <- "B-Sc_2012"
# PRE_ENTRY_1 <- "B-Sc_2012/input/B-Sc_2012_h5 files/PRE_ENTRY_1.h5"
#   hd_file_PRE1 <- PRE_ENTRY_1
# res_PRE_1 <-choosing_sample(hd_file_PRE1, test_orfs, dataset = dataset_Brar, gff_df, min_read_length = 10)
# 
# PRE_ENTRY_2 <-"B-Sc_2012/input/B-Sc_2012_h5 files/PRE_ENTRY_2.h5"
#   hd_file_PRE2 <- PRE_ENTRY_2
# res_PRE_2 <-choosing_sample(hd_file_PRE2, test_orfs, dataset = dataset_Brar, gff_df, min_read_length = 10)
# 
# RECOMB_1 <-"B-Sc_2012/input/B-Sc_2012_h5 files/RECOMB_1.h5"
#   hd_file_RECOMB_1 <- RECOMB_1
# res_RECOMB_1 <- choosing_sample(hd_file_RECOMB_1, test_orfs, dataset = dataset_Brar, gff_df, min_read_length = 10)
#   
# RECOMB_2 <-"B-Sc_2012/input/B-Sc_2012_h5 files/RECOMB_2.h5"
#   hd_file_RECOMB_2 <- RECOMB_2
# res_RECOMB_2 <- choosing_sample(hd_file_RECOMB_2, test_orfs, dataset = dataset_Brar, gff_df, min_read_length = 10)
#   
# SPORES_1 <-"B-Sc_2012/input/B-Sc_2012_h5 files/SPORES_1.h5"
#   hd_file_SPORES_1 <- SPORES_1
# res_SPORES_1 <- choosing_sample(hd_file_SPORES_1, test_orfs, dataset = dataset_Brar, gff_df, min_read_length = 10)
# 
# SPORES_2 <-"B-Sc_2012/input/B-Sc_2012_h5 files/SPORES_2.h5"
#   hd_file_SPORES_2 <- SPORES_2
# res_SPORES_2 <- choosing_sample(hd_file_SPORES_2, test_orfs, dataset = dataset_Brar, gff_df, min_read_length = 10)
# 
# plot_PRE1 <- plotting_5UTR(res_PRE_1)
# plot_PRE2 <- plotting_5UTR(res_PRE_2)
# plot_RECOMB1 <- plotting_5UTR(res_RECOMB_1)
# plot_RECOMB2 <- plotting_5UTR(res_RECOMB_2)
# plot_SPORES1 <- plotting_5UTR(res_SPORES_1)
# plot_SPORES2 <- plotting_5UTR(res_SPORES_2)
# 
# plot_BRAR <- ggarrange(plot_PRE1, plot_PRE2, plot_RECOMB1, plot_RECOMB2, plot_SPORES1, plot_SPORES2, labels = c("PRE-MEIOTIC ENTRY BioRep 1 (0h)","PRE-MEIOTIC ENTRY BioRep 2 (0h)", "RECOMBINATION BioRep 1 (6h)", "RECOMBINATION BioRep 2 (3h)", "     SPORE BioRep 1 (11h)", "     SPORE BioRep 2 (24h)"), ncol = 1, nrow = 6 )
# 


##################################### MISCALLENAEOUS#######################
# ############################################## seems we're not using it at the end but still useful 
# sliding_windows <- function(tibble) {
#   RcppRoll::roll_sum(tibble$Counts, n =30, by = 30) %>% 
#     tibble::enframe(name = NULL) %>%
#      t()
# }
# 
# sliding_windows_multiple <- function(tibble, genes){
#   names(tibble) <- genes
#   joined <-purrr::map(tibble, sliding_windows ) 
# 
#   
#   bind_rows(joined) %>%
#     set_rownames( value = c(seq(from = -250, to = 20, by = 30)))  %>%
#     t() %>%
#     return()  ###add a column instead of the set_rownames and you want these to be ranges (from -250 to -230)
# }
# 
# xxxx <-sliding_windows_multiple(mapped_A_genes, test_orfs)
# 
# selected <- xxxx[xxxx > 10]
# 
# (xxxx > 0) * 1.0 
# 
# if (xxxx > 10 ) {
#   (xxxx*10)
# }
# empty <- list()
# 
# #dziala tylko nie wiem jak wsadzic te values do nowej matrycy wraz ze wspolrzednymi 
# for (i in 1:nrow(xxxx)) {
#    for (j in 2:ncol(xxxx)) {
#     if(xxxx[i,j] >= 10) {
#       
#     }}}
# 
#     # xxxx >0
#     # -250  -220  -190  -160  -130  -100   -70   -40   -10    20
#     # YCR012W FALSE FALSE FALSE FALSE FALSE FALSE  TRUE  TRUE  TRUE  TRUE
#     # YEL009C FALSE FALSE  TRUE  TRUE FALSE FALSE FALSE FALSE FALSE FALSE
#     # YOR303W FALSE FALSE FALSE  TRUE  TRUE  TRUE  TRUE FALSE  TRUE  TRUE
#     # YOL130W FALSE FALSE FALSE FALSE FALSE  TRUE FALSE FALSE FALSE FALSE
#     # YGR094W FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE  TRUE
#     # YML065W FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
# 
# 
# 
# xxxx > 10 
#   
# xxxx_selected <-which( xxxx > 10, arr.ind=T )
# #          row col
# # YEL009C   2   3
# # YOR303W   3   6
# # YCR012W   1   7
# # YOR303W   3   7
# # YCR012W   1   8
# # YCR012W   1   9
# # YCR012W   1  10
# 
# fastQ <- "G-Sc_2014/input/yeast_CDS_w_250utrs.fa"
# c<-read.FASTA(fastQ, type = "DNA")
# # c$YAL068C
# # [1] 88 28 28 18 88 18 48 88 88 88 48 88 18 18 18 88 18 48 88 18 18 28 48 18 18 28 88 48 88 88 88 28 [33] 88 88 48 88 48 28 88 18 28 18 28 28 88 18 88 48 88 48 88 18 88 88 18 48 88 48 88 18 18 48 18 48 [65] 18 48 88 88 88 48 88 18 48 88 48 88 18 88 18 88 48 88 48 88 88 18 88 18 18 18 48 88 18 48 88 48
# 
# b <- read.dna(fastQ, format = "fasta",
#          nlines = 0, comment.char = "#",
#          as.character = TRUE, as.matrix = TRUE)
# # $YDR122W   list
# # [1] "t" "t" "t" "t" "g" "a" "g" "a" "t" "t" "t" "a" "c" "t" "t" "c" "g" "t" "t" "a" "t" "t" "a"
# # [24] "t" "a" "a" "g" "g" "a" "c" "a" "t" "a" "c" "g" "g" "t" "a" "a" "c" "c" "t" "a" "g" "g" "c"
# # [47] "t" "a" "g" "t" "t" "a" "a" "c" "a" "t" "a" "a" "t" "t" "a" "g" "t" "t" "g" "t" "c" "a" "c"
# # [70] "c" "c" "t" "t" "c" "g" "c" "c" "c" "c" "c" "c" "t" "a" "c" "c" "c" "t" "a" "t" "c" "g" "g"
# #