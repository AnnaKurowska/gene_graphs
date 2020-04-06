
############################################################################################
                                        ##set up##
# setwd("/Users/Ania/Desktop/Szkoła/4th year/Dissertation/gene_graphs/")
library(here)
here()
# "/Users/Ania/Desktop/planets" is that the issue for Mrkd?
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
# WT_none <- "G-Sc_2014/output/WTnone/WTnone.h5"
WT_3AT <- "G-Sc_2014/output/WT3AT/WT3AT.h5"
#WT_CHX <- "G-Sc_2014/output/WTCHX/WTCHX.h5"

# prepare files, opens hdf5 file connection
dataset <- "G-Sc_2014"
hd_file <- WT_3AT
hdf5file <- rhdf5::H5Fopen(hd_file) # filehandle for the h5 file

asite_disp_length_file <- "G-Sc_2014/input/asite_disp_length_yeast_standard.txt"

asite_disp_length <- readr::read_tsv(asite_disp_length_file,
                                     comment = "#"
)



#Initial set ups
test_orfs <- c("YCR012W","YEL009C","YOR303W","YOL130W","YGR094W", "YML065W")
gene_names <- rhdf5::h5ls(hdf5file, recursive = 1)$name

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

gff_df <- readGFFAsDf(orf_gff_file) 
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

# Takes value for 3' end of the table (includes whole 5'UTR +50 nt of CDS)
CDS3_end <- function(gene, gffdf) {
  gffdf %>% 
    dplyr::filter(Name == gene, type=="UTR5") %>% 
    dplyr::pull(end) + nnt_gene
} 


# Takes value for length of 5'UTR  
UTR5_length <- function(gene, gffdf) {
  gffdf %>% 
    dplyr::filter(Name == gene, type=="UTR5" ) %>% 
    dplyr::pull(width)
} 

# Creates a matrix with number of columns that start at 5' UTR and finish at 50'th nt of the CDS
GetGeneDatamatrix5UTR <- function(gene, dataset, hdf5file, gffdf, nnt_gene) {
  data_mat_all <- GetGeneDatamatrix(gene, dataset, hdf5file)
  n_left5 <- Get5UTRstart(gene, gffdf) # column to start from (5'end)
  n_right3 <- UTR5_length(gene, gffdf) + nnt_gene # column to end with (3'end) 
  data_mat_5start <- data_mat_all[, n_left5 : n_right3]
  # data_mat_5start <- tibble::as_tibble(data_mat_5start, .name_repair="minimal")
  return(data_mat_5start)
}

# GetPosnCountOutput <- function(x){
#   output_thing <- x %>% 
#     gather(-read_length, key="Position", value = "Counts") %>%
#     group_by(Position) %>%
#     summarise(Counts=sum(Counts)) %>%
#     arrange(as.integer(Position))
#   return(output_thing)
# }

TidyDatamatrix <- function(x, startpos = 1, startlen = 1, gene) {
  # CHECK startpos/off-by-one
  positions <- startpos:(startpos + ncol(x) - 1)
  readlengths <- startlen:(startlen + nrow(x) - 1)
  x %>%
    set_colnames(positions) %>%
    as_tibble() %>%
    mutate(ReadLen = readlengths) %>%
    gather(-ReadLen, key = "Pos", value = "Counts", convert = FALSE) %>%
    mutate(Pos = as.integer(Pos), Counts = as.integer(Counts)) %>%
    group_by(Pos) %>%
    summarise(Counts=sum(Counts))
}

# final function from raw data processing to data visualization
GetGeneDatamatrix5UTR_non_A <- function(gene, dataset, hdf5file, gffdf, nnt_gene) {
  data_mat_all <- GetGeneDatamatrix(gene, dataset, hdf5file)
  n_left5 <- Get5UTRstart(gene, gffdf) # column to start from (5'end)
  n_right3 <- UTR5_length(gene, gffdf) + nnt_gene # column to end with (3'end) 
  data_mat_5start <- data_mat_all[, n_left5 : n_right3]
  data_mat_5start <- tibble::as_tibble(data_mat_5start, .name_repair="minimal")
  return(data_mat_5start)
}

UTR5_table <- function(gene) {
 lapply(gene,
         function(gene)
           GetGeneDatamatrix5UTR_non_A(gene,
                                 dataset,
                                 hdf5file,
                                  gffdf =gff_df,
                                 nnt_gene = 50)
  ) %>%
    Reduce("+", .) %>% # sums the list of data matrices
    TidyDatamatrix(startpos = -250, startlen = 10) %>%
 return()

}

# the final function i'd use :) iteration of UTR5_table
final_function_table <- function(genes){
 table <-purrr::map(genes, UTR5_table)
 names(table) <- genes
 return(table)

}
output_orfs<-final_function_table(test_orfs) #it works just fine 

##meta-analysis:

meta_5genes <- UTR5_table(test_orfs)
  ##just do the plotting here

plotting_meta_analysis<- function(input_data) {
  # text_AUG <- textGrob("AUG", gp=gpar(fontsize=13, fontface="bold")) #to place text annotation

  #the actual plotting
  plotted_UTR <- ggplot(input_data) +
    geom_density(aes(x=Pos, y=Counts), stat="identity") +
    scale_x_continuous(limits = c(-250,50), expand = c(0, 0)) +
    scale_y_continuous(expand = c(0,0)) +
    labs(y= "Read count", x = "Position") +
     ggtitle(paste0("Ribosome footprint density of all genes: WT none")) +
    # coord_cartesian(clip = "off") +
    # annotation_custom(text_AUG,xmin=0,xmax=0,ymin=0,ymax=5) +
    theme_classic() %>%
    return()
}
meta_5genes_plot <- plotting_meta_analysis(meta_5genes)


##I should still be using these functions as they will allow me to compare the A- and non-A site mappng!

#########################A site displacement slot#########################

###
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
  
}  ##I changed n= and by= from 3L to 1L bc I want NT resolution 

# left = Get5UTRstart(test_orfs[1], gff_df) #finds left and right value for each gene 
# right = CDS3_end(test_orfs[1], gff_df)
#potential issue: i used Get5UTRstart and CDS3_end for both GetGeneDataMatrix5UTR and SnapToCodon, which isnt exactly whats done in riboviz i think

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


GetPosition <- function(x, startpos = 1)  {
  positions <- startpos:(startpos + ncol(x) - 1) 
  positions %>% 
    as_tibble() %>%
    set_colnames("Pos")
}  

###
GetGeneCodonPosReads1dsnap <- function(gene, dataset, hdf5file, gffdf , nnt_gene, 
                                       min_read_length, 
                                       asite_disp_length, 
                                       snapdisp = 0L) {
 
  # will have to call x and gff_df the same!
  reads_pos_length <- GetGeneDatamatrix5UTR(gene, dataset, hdf5file,gffdf, nnt_gene)
  reads_asitepos <- CalcAsiteFixed(
    reads_pos_length, min_read_length,
    asite_disp_length)  
  
  Counts <- SnapToCodon(reads_asitepos, left = Get5UTRstart(gene, gff_df) , right = CDS3_end(gene, gff_df), snapdisp)
  Pos <- reads_pos_length %>%
    GetPosition(startpos = -250)
    # TidyDatamatrix(startpos = -250, startlen = 10) %>%
    # select(Pos)

  cbind(Pos, Counts) %>% as_tibble()
}

#  map(test_orfs, ~GetGeneCodonPosReads1dsnap(
#   .,
#   dataset,
#   hdf5file,
#   gff_df,
#   nnt_gene = 50,
#   min_read_length = 10,
#   asite_disp_length = asite_disp_length ))
# names(Counts_Asite_mapped) <- test_orfs

A_mapped_genes <- function(gene, dataset, hdf5file, gffdf, min_read_length) {
  output<- map(gene, GetGeneCodonPosReads1dsnap, dataset, hdf5file, gffdf, nnt_gene, min_read_length, asite_disp_length)
  
  names(output) <- gene
  return(output)
}

A_mapped_genes(test_orfs, dataset = dataset, hdf5file = hdf5file, gffdf = gff_df, min_read_length = 10)

##### META VERSION
  #it would really require the function above to be automated 

# Counts_Asite_mapped  <- map(gene_names, ~GetGeneCodonPosReads1dsnap(
#   .,
#   dataset,
#   hdf5file,
#   gff_df,
#   nnt_gene = 50,
#   min_read_length = 10,
#   asite_disp_length = asite_disp_length ))
# names(Counts_Asite_mapped) <- gene_names
# 
# Counts <-Counts_Asite_mapped %>%
#   Reduce("+", .) %>%
#   select(Counts)
# Pos <- Counts_Asite_mapped$Q0045$Pos
# 
# Counts_Asite_mapped_all <-cbind(Pos, Counts) %>% as_tibble(.name_repair = "minimal")
# 
# plotting_5UTR(Counts_Asite_mapped_all) 
#   %>% return()

All_genesAmapped <- function(gene, dataset, hdf5file, gffdf, min_read_length) {
  output<- map(gene, GetGeneCodonPosReads1dsnap, dataset, hdf5file, gffdf, nnt_gene, min_read_length, asite_disp_length)
  
  names(output) <- gene
  
  Counts <- output %>%
    Reduce("+", .) %>%
    select(Counts)
  Pos <- output$YCR012W$Pos
  
  Counts_Asite_mapped_all <-cbind(Pos, Counts) %>% as_tibble(.name_repair = "minimal") 
  return(Counts_Asite_mapped_all)
  
}

All_genesAmapped(gene = test_orfs, dataset = dataset, hdf5file = hdf5file, gffdf = gff_df, min_read_length = 10)


############################################################################################                                               ##plotting##

plotting_5UTR<- function(input_data) {
  # text_AUG <- textGrob("AUG", gp=gpar(fontsize=13, fontface="bold")) #to place text annotation
  
  #the actual plotting
  plotted_UTR <- ggplot(input_data) +
    geom_density(aes(x=Pos, y=Counts), stat="identity") +
    scale_x_continuous(limits = c(-250,50), expand = c(0, 0)) +
    scale_y_continuous(expand = c(0,0)) +
    labs(y= "Read count", x = "Position") +
    # ggtitle(paste0("Ribosome footprint density of ", input_data)) +
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

####################################################################################
##2 in one##

#this is what i'd like to achieve for the final A-site mapping function 
# choosing_sample <- function(dataset, genes) {
# hd_file <- dataset
# hdf5file <- rhdf5::H5Fopen(hd_file) # filehandle for the h5 file
# 
# final_function_table(genes) %>%
# return()
# }

# choosing_sample(WT_none, genes = test_orfs)

##WT none
WT_none <- "G-Sc_2014/output/WTnone/WTnone.h5"
hd_file <- WT_none
hdf5file <- rhdf5::H5Fopen(hd_file) # filehandle for the h5 file
output_none <-final_function_table(test_orfs)

## WT_3AT
WT_3AT <- "G-Sc_2014/output/WT3AT/WT3AT.h5"
hd_file <- WT_3AT
hdf5file <- rhdf5::H5Fopen(hd_file) # filehandle for the h5 file
output_3AT <-final_function_table(test_orfs)

## WT_CHX
WT_CHX <- "G-Sc_2014/output/WTCHX/WTCHX.h5"
hd_file <- WT_CHX
hdf5file <- rhdf5::H5Fopen(hd_file) # filehandle for the h5 file
output_CHX <-final_function_table(test_orfs)
######

two_in_one_plots <- function(sample1, sample2, gene, cond1_name, cond2_name) {
  ggplot(sample1,aes(x=Pos, y=Counts)) +
    geom_density(data = sample1, stat="identity", aes(color = "WT_3AT")) +
    geom_density(data = sample2, stat="identity", aes(color = "WT_CHX")) +
    scale_x_continuous(limits = c(-100,50), expand = c(0, 0)) +
    scale_y_continuous(expand = c(0,0)) +
    labs(y= "Read count", x = "Position") +
    theme_classic() + 
    ggtitle(paste0(gene, ":", cond1_name, " vs ", cond2_name)) +
    scale_color_manual(values = c(WT_3AT = "red", WT_CHX = "blue"))
}
#### that's the kind of stuff i'd like to achieve 

WT_CHX_3AT <- two_in_one_plots(output_3AT[[1]], output_CHX[[1]], "YCR012W", "WT 3AT", "WT CHX") 
# they're still lists and we want tibbles, otherwise it doesnt work
##I have no idea how to iterate this one

####################################################################################
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
sum_uAUG_multiple<- function(genes, uAUG_start, uAUG_end) {
  #function 
  sapply(genes, sum_uAUG, uAUG_start, uAUG_end) %>%
    return()
  #sapply jest lepsze bo nie tworzy listy ale jakieś problemy są 
    # set_rownames(., paste0(gene_names)) %>%
  #names of columns and rows
}

all_regions_together <- function(genes, gene_names) {
  table_genes <-as.character((paste0(gene_names)))  #why does it make it a factor?
  
  region1 <-sum_uAUG_multiple(genes, uAUG_start = -250, uAUG_end = 0)
  region2 <- sum_uAUG_multiple(genes, uAUG_start = -15, uAUG_end = 10)
  region3 <-sum_uAUG_multiple(genes, uAUG_start = -15, uAUG_end = 0) 
  region4 <-sum_uAUG_multiple(genes, uAUG_start = -60, uAUG_end = 0) 
  region5 <-sum_uAUG_multiple(genes, uAUG_start = -120, uAUG_end = -60)
  
  # region5 <- sum_uAUG_multiple(genes, gene_names, uAUG_start = -15, uAUG_end = 10) %>% unlist(as.numeric())
  
  cbind(table_genes,region1, region2, region3, region4, region5) %>% as_tibble() %>%
    set_colnames(c("genes","-250:0", "-15:10", "-15:-0", "-60:0", "-120:-60")) %>%
    gather( key = "region", value = "count", -genes) %>%
    return()
}

#example: 
together <-all_regions_together(genes = output_orfs,gene_names = test_orfs) 
together$count <- as.numeric(together$count)

Wilcoxon_computed <- function(genes, gene_names) {
  table_genes <-as.character((paste0(gene_names)))  #why does it make it a factor?
  
  region1 <-sum_uAUG_multiple(genes, uAUG_start = -250, uAUG_end = 0)
  region2 <- sum_uAUG_multiple(genes, uAUG_start = -15, uAUG_end = 10)
  region3 <-sum_uAUG_multiple(genes, uAUG_start = -15, uAUG_end = 0) 
  region4 <-sum_uAUG_multiple(genes, uAUG_start = -60, uAUG_end = 0) 
  region5 <-sum_uAUG_multiple(genes, uAUG_start = -120, uAUG_end = -60)
  
  # region5 <- sum_uAUG_multiple(genes, gene_names, uAUG_start = -15, uAUG_end = 10) %>% unlist(as.numeric())
  
  table <-cbind(table_genes,region1, region2, region3, region4, region5) %>% as_tibble() %>%
    set_colnames(c("genes","-250:0", "-15:10", "-15:-0", "-60:0", "-120:-60")) %>%
    gather( key = "genes", value = "count")
  table$count <- as.numeric(table$count)
  
  wilcoxon_computed <- pairwise.wilcox.test(table$count, table$genes, paired = FALSE, p.adjust.method = "none")
  return(wilcoxon_computed)
}
#wilcoxon test -i think its for all genes? although not sure!
Wilcoxon_computed(output_orfs, test_orfs)
# data:  table$count and table$genes 
# 
#          -120:-60 -15:-0 -15:10 -250:0
#   -15:-0 0.69     -      -      -     
#   -15:10 0.69     0.92   -      -     
#   -250:0 0.15     0.35   0.46   -     
#   -60:0  0.69     0.75   1.00   0.35 
# P value adjustment method: none 
# Warning messages:
#   1: In wilcox.test.default(xi, xj, paired = paired, ...) :
#   cannot compute exact p-value with ties
# 2: In wilcox.test.default(xi, xj, paired = paired, ...) :
#   cannot compute exact p-value with ties
# 3: In wilcox.test.default(xi, xj, paired = paired, ...) :
#   cannot compute exact p-value with ties
# 4: In wilcox.test.default(xi, xj, paired = paired, ...) :
#   cannot compute exact p-value with ties


#### pairwise Wilcoxon 

pairwise.wilcox.test(together$count, together$region, p.adjust.method = "none", , alternative = "greater")  

## that gives an overall view of all read counts at each region 
plotted_barplots <- function(data) {
  positions <- c( "-250:0", "-15:10", "-15:-0", "-60:0", "-120:-60")
  
  ggplot(data,aes(x=region, y = count, fill = region)) +
    geom_col() +
    theme_classic() +
    labs(y= "Read count", x = "Region (mRNA)") +
    ggtitle("quantification across all genes") + 
    scale_y_continuous(expand = c(0,0)) +
    scale_x_discrete(limits = positions) 

  ## that gives an overall view of all read counts at each region 
}

## Individual for each gene
plotted_each_barplot <- function(data, gene){
  positions <- c( "-250:0", "-15:10", "-15:-0", "-60:0", "-120:-60")
  value <-filter(data, genes == gene)
  
  ggplot(value,aes(x = region, y = count, fill = region)) +
    geom_col() +
    theme_classic() +
    labs(y= "Read count", x = "Region (mRNA)") +
    ggtitle(paste0(gene)) + 
    scale_y_continuous(expand = c(0,0)) +
    scale_x_discrete(limits = positions)
} 
#WT none
together <-all_regions_together(genes = output_none,gene_names = test_orfs) 
together$count <- as.numeric(together$count)

plotted_each_barplot(together, "YCR012W")
plotted_barplots(together)

#WT 3AT
together_3AT <-all_regions_together(genes = output_3AT,gene_names = test_orfs) 
together_3AT$count <- as.numeric(together_3AT$count)

plotted_each_barplot(together_3AT, "YCR012W")
plotted_barplots(together_3AT)

#WT CHX
together_CHX <-all_regions_together(genes = output_CHX,gene_names = test_orfs) 
together_CHX$count <- as.numeric(together_CHX$count)

plotted_each_barplot(together_CHX, "YCR012W")
plotted_barplots(together_CHX)

  # I don't know what to do about it, it's clearly an issue with subsetting 

#example 
apply(together, 1, plotted_each_barplot, test_orfs)
## Error in UseMethod("filter_") : 
##no applicable method for 'filter_' applied to an object of class "list"

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

##########################################################################################
                                     ######FASTQ######
fq.file <- "/Users/Ania/Desktop/Szkoła/4th year/Dissertation/gene_graphs/Ania_Riboviz/G-Sc_2014/input/yeast_CDS_w_250utrs.fa"

fastq_sequence <- Fastqfile(fq.file)
  #not entirely what i want, a bit confusing: 
# A tibble: 37,762 x 3
# Header                        Sequence                      Quality                   
# <chr>                         <chr>                         <chr>                       
#   1 >YAL068C                      ACCTATGAAAGATTTATGATTCGTTCAG… GCTGCTTCAACTATATGCCTTTGAGAAT…
# 2 AACAAATACAATGGTCAAATTAACTTCA… TAGCTCAATCTGACGAAAGAGTCAACTT… CACCATGTTGACCGGTATTGCTCCAGAC…
# 3 CAGCCATCTCCAGTGCTCTATCCAAGGA… TTCCATAGAAATTGAAAATTAACGAACA… AAAGAAACTTCTACACTATTGTAGAAAA…
# 4 >YAL067W-A                    ATATTCTCAAAGGCATATAGTTGAAGCA… TTCTGAACGAATCATAAATCTTTCATAG…
# 5 GAATGTGGGAATGCCAATTATAGGGGTG… TCAAAAAGAATATCCGAATTTTAGATTT… TCTGTAAACTTGTGAACTCTCGGCAAAT…
# 6 CAAGTTGATATCAAACAGATACATATTT… ACGGCTAACTGAACCTAAGTAGGGATAT… >YAL067C                     
# 7 AAAAAAACTGCTACAAATAATAAATAAA… ACCAAGCGGAGACATGCGTTTAGATGAG… AATAACATACATGTATTCAATTGTTAAA…
# 8 AGCGGCAGGTGGAAGACCTGCCAGATGA… ACAGCTGAAAATTTCATCACGACTACAA… ATCGGATCAATGAAAAGGAAAGATCTCA…
# 9 TTATTAATTAAATTGGATGTCCTTTTAG… AAACAACGCTTACGTTTCGGGAATGAAG… GACTTATGTTGGTCGCTTTTAACCGTTG…
# 10 TGGGGCTTTTGAAGCGCCAAGTTATTTG… GTTCTGCTTTTTACTATTTGGGCCAGTA… ATTTTACTCCCTGCCAGGTGACCCATAC…
# … with 37,752 more rows
