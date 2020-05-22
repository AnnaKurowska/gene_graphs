
# set up the directory 
library(here)
here()

# load the packages
library(rhdf5)
library(tidyr)
library(dplyr)
library(magrittr)
library(purrr)
library(ggplot2)
library(ggpubr)
library(tibble)
library(grid)
library(ggrepel)


################################### Functions ###################################

#Creates a gff table  
readGFFAsDf <- purrr::compose(
  rtracklayer::readGFFAsGRanges,
  data.frame, 
  as_tibble,
  .dir = "forward" # functions called from left to right
)

# Get the total number of reads from each gene 
GetGeneReadsTotal <- function(gene, dataset, hdf5file) {
  rhdf5::H5Aread(rhdf5::H5Aopen(rhdf5::H5Gopen(hdf5file, paste0("/", gene, "/", dataset, "/reads")), "reads_total"))
}

#Creates a matrix from hdf5file
GetGeneDatamatrix <- function(gene, dataset, hdf5file) {
  hdf5file %>%
    rhdf5::H5Dopen(
      name = paste0("/", gene, "/", dataset, "/reads/data")
    ) %>%
    rhdf5::H5Dread() %>%
    return()
}

# Takes value for start position of TL for matrix creation (no negative values e.g.: 1)
GetTLstart <- function(gene, gffdf) {
  gffdf %>% 
    dplyr::filter(Name == gene, type=="UTR5" ) %>% 
    dplyr::pull(start)
}

# Takes value for 3' end of the table (includes whole TL +50 nt of CDS)
CDS3_end <- function(gene, gffdf) {
  gffdf %>% 
    dplyr::filter(Name == gene, type=="UTR5") %>% 
    dplyr::pull(end) + nnt_gene
} 

# Creates a matrix with number of columns that start at 5' TL and finish at 50'th NT of  CDS
GetGeneDatamatrixTL <- function(gene, dataset, hdf5file, gffdf, nnt_gene) {
  data_mat_all <- GetGeneDatamatrix(gene, dataset, hdf5file)
  n_left5 <- GetTLstart(gene, gffdf) # column to start from (5'end)
  n_right3 <-CDS3_end(gene, gffdf)  # column to end with (3'end)
  data_mat_5start <- data_mat_all[, n_left5 : n_right3]
  # data_mat_5start <- tibble::as_tibble(data_mat_5start, .name_repair="minimal")
  return(data_mat_5start)
}
#### functions for A-site mapping
## A-site mapping function. Calculate read A-site using a fixed displacement for a single read length
CalcAsiteFixedOneLength <- function(reads_pos_length, min_read_length,
                                    read_length, asite_disp) {
  
  length_row_choose <- read_length - min_read_length + 1
  reads_pos_length[length_row_choose, ] %>%
    dplyr::lag(n = asite_disp, default = 0)
}


## A-site mapping function. Calculate read A-site using a fixed displacement for fixed read lengths
CalcAsiteFixed <- function(reads_pos_length, min_read_length,
                           asite_disp_length,
                           colsum_out = TRUE) {
  
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

GetPosition <- function(x, startpos = 1)  {
  positions <- startpos:(startpos + ncol(x) - 1) 
  positions %>% 
    tibble::as_tibble() %>%
    magrittr::set_colnames("Pos")
}  

## Final A-site mapping function
GetGene <- function(gene, dataset, hdf5file, gffdf,
                    nnt_gene, min_read_length, asite_disp_length,  snapdisp = 0L, scaling_factor) {
  
  left  = GetTLstart(gene, gffdf)
  right = CDS3_end(gene, gffdf) 

  reads_pos_length <- GetGeneDatamatrixTL(gene, dataset,
                                            hdf5file, gffdf,nnt_gene)
  
  
  Counts <- CalcAsiteFixed(reads_pos_length,
                                   min_read_length,
                                   asite_disp_length)  
  
  
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

##A-site mapping for multiple genes; output is an individual tibble for each gene
GetGeneMultiple <- function(gene,
                            dataset,
                            hdf5file,
                            gffdf,
                            min_read_length, 
                            asite_disp_length,
                            scaling_factor) {
  
  
  output<- purrr::map(gene,
                      GetGene,
                      dataset,
                      hdf5file,
                      gffdf,
                      nnt_gene,
                      min_read_length,
                      asite_disp_length,
                      scaling_factor = 2.418605)
  
  names(output) <- gene
  
  
  return(output)
}


## A-site mapping for multiple genes; output is one tibble with cumulative count of reads at each position
GetGeneMeta <- function(gene, dataset, hdf5file, gffdf,
                        min_read_length, asite_disp_length,
                        scaling_factor) {
  
  output<- purrr::map(gene,
                      GetGene,
                      dataset,
                      hdf5file,
                      gffdf,
                      nnt_gene,
                      min_read_length,
                      asite_disp_length,
                      scaling_factor = 2.41)
  
  names(output) <- gene
  
  Counts <- output %>%
    Reduce("+", .) %>%
    dplyr::select(Counts)
  
  Pos <- output$YEL009C$Pos 
  
  Counts_Asite_mapped_all <-base::cbind(Pos, Counts) %>%
    tibble::as_tibble(.name_repair = "minimal")
  
  return(Counts_Asite_mapped_all)
  
}

## Combines 2 tables from GetGene or GetGeneMeta into one table 
join_GetGeneData <-function(dataset1, dataset2){
  full_join(dataset1,dataset2, by = "Pos") %>%
    mutate(., Counts = Counts.x + Counts.y) %>%
    select(Pos, Counts)%>%
    return()
} 

#### functions which take output from GetGene and GetGeneMultiple as input

# Density in TL (-250:-5) and AUG (-5: 10)
CanVsNon <- function(dataset) {
  
  TL <- bind_rows(dataset, .id = "Gene") %>%
    dplyr::filter( Pos >= -250, Pos < -5) %>%
    dplyr::group_by(Gene) %>%
    dplyr::summarise(Count_TL = sum(Counts)/(length(-250:-5)-1))
  
  AUG <- bind_rows(dataset, .id = "Gene") %>%
    dplyr::filter( Pos >= -5, Pos < 10) %>%
    dplyr::group_by(Gene) %>%
    dplyr::summarise(Count_AUG = sum(Counts)/(length(-5:10)-1))
  
  TL_AUG <-dplyr::full_join(TL, AUG, by = "Gene") 
  return(TL_AUG)
  
}

###Scatter plot showing the TL/AUG ratio, input: CanVsNon output
TLvsAUGscatter <-function(CanNonData){
  ggplot(CanNonData, aes(x = Count_AUG, y = Count_TL)) +
    geom_point() +
    scale_y_log10() +
    scale_x_log10() +
    geom_text_repel(aes(label = Gene)) +
    labs(x= "AUG (CPMB)", y = "TL (CPMB)") +
    theme_classic()
  
}

##### TL density for same gene across different stages, input: output of CanVsNon data for each stage

TL_AUG_Regulation <-function(CanNonData1, CanNonData2, CanNonData3){
  
  comparing_TL_test <- full_join(CanNonData1, CanNonData2, by = "Gene") %>%
    full_join(CanNonData3, by = "Gene") %>%
    set_colnames(c("Gene", "TL_WT_NONE", "AUG_WT_NONE", "TL_WT_CHX","AUG_WT_CHX", "TL_WT_3AT", "AUG_WT_3AT" )) 
  
  comparing_TL_test <-comparing_TL_test%>%
    select("Gene","TL_WT_NONE","TL_WT_CHX","TL_WT_3AT") %>%
    gather(key = "Condition", value = "Reads", -Gene)
  
  comparing_TL_test$Condition <- factor(comparing_TL_test$Condition, level=unique(comparing_TL_test$Condition))  
  
  ggplot(comparing_TL_test, aes(x = Condition, group = Gene, y= Reads, color = Gene)) +
    geom_point() +
    geom_line() +
    scale_y_continuous(expand = c(0,0)) +
    geom_text_repel(aes(label = round(Reads, digits = 3)))
  
}

## distribution plot, input: output from GetGene or GetGeneMeta

DistributionPlot<- function(input_data) {
  
  plotted_TL <- ggplot2::ggplot(input_data) +
    geom_density(aes(x=Pos, y=Counts), stat="identity") +
    scale_y_continuous(expand = c(0,0)) +
    scale_x_continuous(breaks = c(-250,-200,-150,-100,-60,-30,0,30,50), label=c("-250","200","-150","100", "-60","-30","0\nAUG",30,50)) +
    labs( x = "Position", y = "CPM") +
    theme_classic() 
  return(plotted_TL)
}


## Multiple distribution plots, input: output from GetGeneMultiple, GetGene or GetGeneMeta

MultipleDistributionPlots <- function(input_data, gene_names) {
  purrr :: map(input_data, DistributionPlot) %>%
    ggarrange(plotlist = ., label = gene_names) %>%
    return()
}


## distribution plot of 2 samples, input: GetGene,GetGeneMultiple (1 gene only) or GetGeneMeta

TwoPlotsDistribution <- function(sample1, sample2, gene, cond1_name, cond2_name) {
  ggplot2::ggplot(sample1,aes(x=Pos, y=Counts, alpha = 0.5) ) +
    geom_density(data = sample1, stat="identity", aes(color = cond1_name, fill = cond1_name)) +
    geom_density(data = sample2, stat="identity", aes(color = cond2_name, fill = cond2_name)) +
    scale_x_continuous(limits = c(-250,50), expand = c(0, 0)) +
    scale_y_continuous(expand = c(0,0)) +
    labs(y= "CPM", x = "Position") +
    theme_classic() +
    ggtitle(paste0(gene, ":", cond1_name, " vs ", cond2_name)) +
    theme(legend.title = element_blank())
}

## distribution plot of 3 samples, input: GetGene,GetGeneMultiple (1 gene only) or GetGeneMeta

ThreePlotsDistribution <- function(sample1, sample2, sample3, cond1_name, cond2_name, cond3_name, gene) {
  ggplot2::ggplot(sample1,aes(x=Pos, y=Counts, alpha = 0.5)) +
    geom_density(data = sample1, stat="identity", aes(color = cond1_name, fill = cond1_name )) +
    geom_density(data = sample2, stat="identity", aes(color = cond2_name, fill = cond2_name)) +
    geom_density(data = sample3, stat="identity", aes(color = cond3_name, fill =cond3_name )) +
    scale_x_continuous(breaks = c(-250,-200,-150,-100,-60,-30,0,30,50), label=c("-250","200","-150","100", "-60","-30","0\nAUG",30,50)) +
    scale_y_continuous(expand = c(0,0)) +
    labs(y= "CPM", x = "Position") +
    theme_classic() +
    scale_color_manual(values = c(VEG = "red", ANAPH_I = "blue", SPORE = "black" ))
}


##Sums read counts and normalizes them for density at a selected region for 1 gene, input GetGene, GetGeneMeta
sum_uAUG <- function(datatibble, uAUG_start, uAUG_end){
  datatibble %>%
    dplyr::filter(Pos >= uAUG_start, Pos <= uAUG_end) %>%
    dplyr::summarise(sum_read_counts_uAUG=sum(Counts)/(length(uAUG_start:uAUG_end)-1)) 
}

##Sums read counts and normalizes them for density at a selected region for many genes, input GetGene, GetGeneMultiple, GetGeneMeta
sum_uAUG_multiple<- function(genes, uAUG_start, uAUG_end) {
  sapply(genes, sum_uAUG, uAUG_start, uAUG_end) %>%
    return()
}


### Used in GetTLvsAUG. Calculates footprint density in the TL and AUG, input GetGene, GetGeneMultiple or GetGeneMeta. 

CDS_TL <- function(genes, gene_names) {
  table_genes <-as.character((paste0(gene_names)))  #why does it make it a factor?
  
  
  region1 <-sum_uAUG_multiple(genes, uAUG_start = -250, uAUG_end = -5)
  region2 <- sum_uAUG_multiple(genes, uAUG_start = -5, uAUG_end = 10)
  
  
  result <-base::cbind(table_genes, region1, region2) %>%
    as_tibble() %>%
    set_colnames(c("genes","TL", "AUG"))
  
  result_new <- result %>%
    tidyr::gather( key = "region", value = "count", -genes)
  result_new$count <- as.numeric(result_new$count)
  return(result_new)
}


# efficiency bar plot (TL and AUG) for all input genes, includes the CDS_TL step; input: GetGene, GetGeneMultiple or GetGeneMeta
GetTLvsAUG <- function(genes, gene_names) {
  CDS_TL_data <-CDS_TL(genes, gene_names) %>%
    tidyr::gather( key = "region", value = "count", -genes)
  CDS_TL_data$count <- as.numeric(CDS_TL_data$count)
  
  positions <- c("TL", "AUG")
  
   barplots<-ggplot(CDS_TL_data, aes(x=region, y = count, fill = region )) +
    geom_col(show.legend = FALSE) +
    theme_classic() +
    labs(y= "CPMB", x = "Region (mRNA)") +
    scale_y_continuous(expand = c(0,0)) +
    scale_x_discrete(limits = positions) 
  return(barplots)
}

GetTLvsAUG(genes = exemplary_genes_anaphase[[1]],gene_names = exemplary_genes[[1]])

#Creates a table with read density at different regions, input:GetGene,GetGeneMultiple (1 gene only) or GetGeneMeta

UpstreamRegionsCount <- function(genes, gene_names) {
  table_genes <-as.character((paste0(gene_names)))  #why does it make it a factor?
  
  region1 <-sum_uAUG_multiple(genes, uAUG_start = -250, uAUG_end = 0)
  region2 <- sum_uAUG_multiple(genes, uAUG_start = -5, uAUG_end = 10)
  region3 <-sum_uAUG_multiple(genes, uAUG_start = -15, uAUG_end = 0)
  region4 <-sum_uAUG_multiple(genes, uAUG_start = -60, uAUG_end = 0)
  region5 <-sum_uAUG_multiple(genes, uAUG_start = -120, uAUG_end = -60)
  region6  <-sum_uAUG_multiple(genes, uAUG_start = -250, uAUG_end = -120)
  
  
  TL_regions <-base::cbind(table_genes, region1, region2, region3, region4, region5, region6) %>%
    as_tibble() %>%
    magrittr::set_colnames(c("genes","TL", "AUG", "-15:-0", "-60:0", "-120:-60", "-250:-60")) %>%
    tidyr::gather( key = "region", value = "read", -genes)
  TL_regions$read <- as.numeric(TL_regions$read)
  
  return(TL_regions)
}


##Bar plots showing read density at each region for all selected genes (all or multiple), input: output of all_regions_together
Different_Regions_barplot <- function(data) {
  
  ggplot(data, aes(x = region, y = read, fill = region)) +
    geom_col(show.legend = FALSE) +
    theme_classic() +
    labs(y= "CPMB", x = "Region (mRNA)") +
    ggtitle("quantification across all genes") +
    scale_y_continuous(expand = c(0,0))
  # scale_x_discrete(limits = positions)
  
}

##Bar plots showing read density at each region for one gene, input: output of all_regions_together

plotted_each_barplot <- function(data, gene) {
  value <-filter(data, genes == gene)
  
  ggplot(value, aes(x = region, y = read, fill = region)) +
    geom_col(show.legend = FALSE) +
    theme_classic() +
    labs(y= "CPMB", x = "Region (mRNA)") +
    ggtitle(paste0(gene)) +
    scale_y_continuous(expand = c(0,0)) 
  # scale_x_discrete(limits = positions)
}

## Comparison of TL to AUG for all or multiple selected genes in 1 condition, input: Output of CanVsNon

all_genes_TL_AUG <-function(condition1_TL){
  TL <- sum(condition1_TL$Count_TL)
  AUG <- sum(condition1_TL$Count_AUG)
  
  Condition <- c("TL", "AUG")
  for_plot_TL <- rbind(TL, AUG) %>%
    as_tibble %>%
    set_colnames("TL") %>%
    cbind(Condition)
  
  plot_TL <-ggplot(for_plot_TL) +
    geom_bar(aes(x=Condition, y= TL, fill = Condition), stat="identity", show.legend = FALSE ) +
    scale_y_continuous(expand = c(0,0)) +
    theme_classic() +
    labs(x= "Region", y = "TL (CPMB) ")
  return(plot_TL)
}

### Comparison of TL density for all genes between 3 conditions, input: output of CanVsNon 

TL_3_conditions_all <-function(condition1_TL, condition2_TL, condition3_TL, name1, name2, name3){
  condition1 <- sum(condition1_TL$Count_TL)
  condition2 <- sum(condition2_TL$Count_TL)
  condition3 <- sum(condition3_TL$Count_TL)
  
  Condition <- paste0(c(name1, name2, name3))
  # Condition <- c("WT NONE", "WT CHX", "WT 3AT")
  for_plot_TL <- rbind(condition1, condition2, condition3) %>%
    as_tibble %>%
    set_colnames("TL") %>%
    cbind(Condition)
  
  for_plot_TL$Condition <- factor(for_plot_TL$Condition, levels=unique(for_plot_TL$Condition))  
  
  plot_TL <-ggplot(for_plot_TL) +
    geom_bar(aes(x=Condition, y= TL, fill = Condition), stat="identity", show.legend = FALSE ) +
    # scale_x_discrete(labels=c("WT_3AT","WT_CHX", "WT_NONE")) +
    scale_y_continuous(expand = c(0,0)) +
    theme_classic() +
    labs(x= "Stage", y = "TL (CPMB)")
  return(plot_TL)
}


TL_3_conditions_all(gene1_scatter,gene2_scatter, gene3_scatter, name1 = "ves", "x", "z")

### Barplot Comparison of TL for 1 gene between 3 conditions, input: output of CanVsNon  
compared_TL_single<-function(condition1, condition2, condition3, gene){
  
  search_for_TL <- function(condition, gene){
    condition %>%
      filter(Gene == gene ) %>%
      select(Count_TL) %>%
      pull
  }
  
  TL1 <-search_for_TL(condition1, gene)
  
  TL2 <-search_for_TL(condition2, gene)
  
  TL3 <-search_for_TL(condition3, gene)
  
  Condition <- c("WT NONE", "WT CHX", "WT 3AT")

  for_plot <- rbind(TL1,TL2,TL3) %>%
    # set_colnames("TL") %>%
    cbind(Condition) %>%
    as_tibble() %>%
    set_colnames(c("x", "TL", "Condition")) %>%
    select(TL, Condition) 
  for_plot$TL <-as.numeric(for_plot$TL)

  plot_TL <-ggplot(for_plot) +
    geom_bar(aes(x=Condition, y= TL, fill = Condition), stat="identity",show.legend = FALSE ) +
    scale_y_continuous(expand = c(0,0)) +
    theme_classic() +
    labs(x= "Condition", y = "CPMB TL") +
    ggtitle(paste0(gene))

  return(plot_TL)
  
}

##function used within TL_AUG_plot_correlation and TL_AUG_correlation; input: output of CanVsNon

Efficiency_preparation<- function(condition1, condition2, condition3){
  
  comparing_genes_stages <- full_join(condition1, condition2, by = "Gene") %>%
    full_join(condition3, by = "Gene") %>%
    set_colnames(c("Gene", "VEG_TL", "VEG_AUG", "ANAPH_I_TL","ANAPH_I_AUG", "SPORE_TL","SPORE_AUG" )) 
  
  comparing_genes_stages_TL <-comparing_genes_stages %>%
    select("Gene","VEG_TL","ANAPH_I_TL","SPORE_TL") %>%
    set_colnames(c("Gene", "VEG", "ANAPH_I", "SPORE")) %>%
    gather(key = "Condition", value = "Reads_TL", -Gene) 
  
  comparing_genes_stages_AUG <-comparing_genes_stages %>%
    select("Gene","VEG_AUG","ANAPH_I_AUG","SPORE_AUG") %>%
    set_colnames(c("Gene", "VEG", "ANAPH_I", "SPORE")) %>%
    gather(key = "Condition", value = "Reads_AUG", -Gene)
  
  combined_TL_AUG <- right_join(comparing_genes_stages_TL, 
                                comparing_genes_stages_AUG, by = c("Gene","Condition"))
  
  return(combined_TL_AUG)
}
somethin2_scatter
### Outputs a figure with TL and AUG footprint densities across 3 conditions, input: output of CanVsNon. used for 1 gene only

Efficiency_AUG_TL <- function(condition1, condition2, condition3, gene, condition_name1,condition_name2, condition_name3 ){
  
  C1 <-filter(condition1, Gene == gene)
  C2 <-filter(condition2, Gene == gene)
  C3 <-filter(condition3, Gene == gene)
  
  prepared_table <- Efficiency_preparation(C1, C2, C3)
  
  prepared_table$Condition <- factor(prepared_table$Condition, levels=unique(prepared_table$Condition))
  
  plot <-ggplot(prepared_table) +
    geom_point(aes(x = Condition, group = Gene, y= Reads_AUG, color = "AUG" )) +
    geom_line(aes(x = Condition, group = Gene, y= Reads_AUG, color = "AUG")) +
    geom_point(aes(x = Condition, group = Gene, y= Reads_TL, color = "TL")) +
    geom_line(aes(x = Condition, group = Gene, y= Reads_TL, color = "TL")) +
    scale_y_log10() +
    labs(y = "CPMB", x = "Condition") +
    guides(color=guide_legend(NULL)) +
    scale_x_discrete(labels= c(paste0(condition_name1), paste0(condition_name2), paste0(condition_name3))) +
    theme_classic()
  return(plot)
  
}

##calculates correlation between the TL and AUG densities across 3 conditions, input: output of CanVsNon
TL_AUG_correlation <- function(condition1, condition2, condition3){
  
  prepared_table <- Efficiency_preparation(condition1, condition2, condition3)
  
  cor.test(prepared_table$Reads_AUG, prepared_table$Reads_TL, method=c("pearson")) %>%
    return()
}

# ################### A SITE MAPPING #########

## to be used in GetGeneDatamatrixTL_non_A
TidyDatamatrix <- function(x, startpos = 1, startlen = 1, gene) {
  # CHECK startpos/off-by-one
  positions <- startpos:(startpos + ncol(x) - 1)
  readlengths <- startlen:(startlen + nrow(x) - 1)
  x %>%
    magrittr::set_colnames(positions) %>%
    tibble::as_tibble() %>%
    dplyr::mutate(ReadLen = readlengths) %>%
    tidyr::gather(-ReadLen, key = "Pos", value = "Counts", convert = FALSE) %>%
    dplyr::mutate(Pos = as.integer(Pos), Counts = as.integer(Counts)) %>%
    dplyr::group_by(Pos) %>%
    dplyr::summarise(Counts=sum(Counts))
}
# 
# final function from raw data processing to data visualization
GetGeneDatamatrixTL_non_A <- function(gene, dataset, hdf5file, gffdf, nnt_gene) {
  data_mat_all <- GetGeneDatamatrix(gene, dataset, hdf5file)
  n_left5 <- GetTLstart(gene, gffdf) # column to start from (5'end)
  n_right3 <- TL_length(gene, gffdf) + nnt_gene # column to end with (3'end)
  data_mat_5start <- data_mat_all[, n_left5 : n_right3]
  data_mat_5start <- tibble::as_tibble(data_mat_5start, .name_repair="minimal")
  return(data_mat_5start)
}

TL_table <- function(gene) {
  lapply(gene,
         function(gene)
           GetGeneDatamatrixTL_non_A(gene,
                                       dataset,
                                       hdf5file,
                                       gffdf =gff_df,
                                       nnt_gene = 50)
  ) %>%
    Reduce("+", .) %>% # sums the list of data matrices
    TidyDatamatrix(startpos = -250, startlen = 10) %>%
    return()
  
}


## iteration of TL_table over multiple genes
final_function_table <- function(genes){
  table <-purrr::map(genes, TL_table)
  names(table) <- genes
  return(table)
  
}

plotting_meta_analysis<- function(input_data) {
  
  plotted_TL <- ggplot(input_data) +
    geom_density(aes(x=Pos, y=Counts), stat="identity") +
    scale_x_continuous(limits = c(-250,50), expand = c(0, 0)) +
    scale_y_continuous(expand = c(0,0)) +
    labs(y= "Read count", x = "Position") +
    ggtitle(paste0("Ribosome footprint density of all genes: WT 3-AT")) +
    theme_classic() %>%
    return()
}

################################### Set up ######################################


#set up for Guydosh and Green (2014)
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

# set up for Brar et al. (2012) 
#VEG_1

VEG_1 <- "B-Sc_2012/input/B-Sc_2020_master_20200428_6samples/output/VEG_1/VEG_1.h5"
hd_file_VEG_1 <- VEG_1
hdf5file_VEG_1 <- rhdf5::H5Fopen(hd_file_VEG_1)

#VEG_2

VEG_2 <- "B-Sc_2012/input/B-Sc_2020_master_20200428_6samples/output/VEG_2/VEG_2.h5"
hd_file_VEG_2 <- VEG_2
hdf5file_VEG_2 <- rhdf5::H5Fopen(hd_file_VEG_2)

# ANAPH_I_1
ANAPH_1 <- "B-Sc_2012/input/B-Sc_2020_master_20200428_6samples/output/ANAPH_1/ANAPH_1.h5"
hd_file_A_1 <- ANAPH_1
hdf5file_A_1 <- rhdf5::H5Fopen(hd_file_A_1)

# ANAPH_I_2
ANAPH_2 <- "B-Sc_2012/input/B-Sc_2020_master_20200428_6samples/output/ANAPH_2/ANAPH_2.h5"
hd_file_A_2 <- ANAPH_2
hdf5file_A_2 <- rhdf5::H5Fopen(hd_file_A_2)

#SPORE_1

SPORE_1 <- "B-Sc_2012/input/B-Sc_2020_master_20200428_6samples/output/SPORE_1/SPORE_1.h5"
hd_file_SPORE_1<- SPORE_1
hdf5file_SPORE_1 <- rhdf5::H5Fopen(hd_file_SPORE_1)

#SPORE_2

SPORE_2 <- "B-Sc_2012/input/B-Sc_2020_master_20200428_6samples/output/SPORE_2/SPORE_2.h5"
hd_file_SPORE_2<- SPORE_2
hdf5file_SPORE_2 <- rhdf5::H5Fopen(hd_file_SPORE_2)
########

asite_disp_length_file <- "G-Sc_2014/input/asite_disp_length_yeast_standard.txt"

asite_disp_length <- readr::read_tsv(asite_disp_length_file,
                                     comment = "#"
)

dataset_G2014 <- "G-Sc_2014"
dataset_Brar <- "B-Sc_2012"
genes <- c("YPR036W-A", "YEL009C", "YOR303W", "YGR094W", "YOL130W", "YOR061W", "YAL040C", "YKL109W", "YIL144W", "YGR037C", "YJL106W", "YBR160W", "YDR172W", "YML065W", "YBR257W", "YOL104C", "YJR094C","YOL125W")
gene_names <- rhdf5::h5ls(hdf5file_none, recursive = 1)$name
exemplary_genes <- c("YCR012W","YEL009C","YOR303W","YOL130W","YGR094W", "YML065W")

nnt_gene<- 50
startpos <-250
startlen <- 10
min_read_count <- 10
orf_gff_file <- "G-Sc_2014/input/yeast_CDS_w_250utrs.gff3"

## Finding the scaling factor:

scaling_factor_none <<- sapply(gene_names, GetGeneReadsTotal, dataset_G2014, hdf5file_none) %>% sum()/ 10^6

scaling_factor_CHX <<- sapply(gene_names, GetGeneReadsTotal, dataset_G2014, hdf5file_CHX) %>% sum()/ 10^6

scaling_factor_3AT <<- sapply(gene_names, GetGeneReadsTotal, dataset_G2014, hdf5file_3AT) %>% sum()/ 10^6

####

scaling_factor_VEG_1 <<- sapply(gene_names, GetGeneReadsTotal, dataset_Brar, hdf5file_VEG_1) %>% sum()/ 10^6

scaling_factor_VEG_2 <<- sapply(gene_names, GetGeneReadsTotal, dataset_Brar, hdf5file_VEG_2) %>% sum()/ 10^6

##

scaling_factor_A_1 <<- sapply(gene_names, GetGeneReadsTotal, dataset_Brar, hdf5file_A_1) %>% sum()/ 10^6

scaling_factor_A_2 <<- sapply(gene_names, GetGeneReadsTotal, dataset_Brar, hdf5file_A_2) %>% sum()/ 10^6


##

scaling_factor_SPORE_1 <<- sapply(gene_names, GetGeneReadsTotal, dataset_Brar, hdf5file_SPORE_1) %>% sum()/ 10^6

scaling_factor_SPORE_2 <<- sapply(gene_names, GetGeneReadsTotal, dataset_Brar, hdf5file_SPORE_2) %>% sum()/ 10^6

####################### Example for each function
# readGFFAsDF:
gff_df <- readGFFAsDf(orf_gff_file) 
# seqnames  start   end width strand source   type  score phase Name     
# <fct>     <int>  <int> <int> <fct>  <fct>    <fct> <dbl> <int> <chr>    
# 1 YAL068C   1   250   250 +      rtracklayer UTR5   NA    NA YAL068C  
# 2 YAL068C   251 613   363 +      rtracklayer CDS    NA    NA YAL068C 

#GetGeneDatamatrix
GetGeneDatamatrix(gene = exemplary_genes[1], dataset = dataset_Brar,hdf5file = hdf5file_SPORE_1)
# [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12] [,13]
# [,14] [,15] [,16] [,17] [,18] [,19] [,20] [,21] [,22] [,23] [,24] [,25]

#GetTLstart
GetTLstart(gene = exemplary_genes[1],gffdf = gff_df)
# [1] 1

#CDS3_end
CDS3_end(gene = exemplary_genes[1],gffdf = gff_df)
# [1] 300

#GetGeneDatamatrixTL
GetGeneDatamatrixTL(gene =exemplary_genes[1],dataset = dataset_Brar,hdf5file = hdf5file_SPORE_1,gffdf = gff_df,nnt_gene = nnt_gene )
# [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12] [,13]
# [1,]    0    0    0    0    0    0    0    0    0     0     0     0     0
# [2,]    0    0    0    0    0    0    0    0    0     0     0     0     0
# [3,]    0    0    0    0    0    0    0    0    0     0     0     0     0

#CalcAsiteFixedOneLength
CalcAsiteFixedOneLength(reads_pos_length = GetGeneDatamatrix(gene = exemplary_genes[1],dataset = dataset_Brar,hdf5file = hdf5file_SPORE_1), min_read_length = 10, read_length = 28, asite_disp = 15)
# [937] 0 0 0 0 0 0 0 0 0 0 3 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
# [973] 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
# [ reached getOption("max.print") -- omitted 751 entries ]

#CalcAsiteFixed
CalcAsiteFixed(reads_pos_length = GetGeneDatamatrix(gene = exemplary_genes[1], dataset = dataset_Brar, hdf5file = hdf5file_SPORE_1), asite_disp_length = asite_disp_length, min_read_length = 10, colsum_out = TRUE)
# [937] 0 0 0 0 0 0 0 3 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
# [973] 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
# [ reached getOption("max.print") -- omitted 751 entries ]

#GetPosition
GetPosition(x =GetGeneDatamatrixTL(exemplary_genes[1], dataset_Brar, hdf5file_SPORE_1, gff_df, nnt_gene))
# Pos
# <int>
#   1     1
# 2     2
# 3     3
# 4     4
# 5     5
# 6     6
# 7     7
# 8     8
# 9     9
# 10    10
# # … with 290 more rows

#GetGene
GetGene(gene = exemplary_genes[1], dataset = dataset_Brar, hdf5file = hdf5file_SPORE_1,gffdf = gff_df,nnt_gene = nnt_gene,min_read_length = 10,asite_disp_length = asite_disp_length,snapdisp = 0L,scaling_factor = scaling_factor_SPORE_1)
# A tibble: 300 x 2
#     Pos Counts
#     <int>  <dbl>
# 1  -250      0
# 2  -249      0
# 3  -248      0
# 4  -247      0
# 5  -246      0
# 6  -245      0
# 7  -244      0
# 8  -243      0
# 9  -242      0
# 10  -241      0
# # … with 290 more rows

#GetGeneMultiple
exemplary_genes_vegetative <-GetGeneMultiple(gene = exemplary_genes,dataset = dataset_Brar,hdf5file = hdf5file_VEG_1,gffdf = gff_df, min_read_length = 10, asite_disp_length = asite_disp_length,, scaling_factor = scaling_factor_VEG_2)

exemplary_genes_anaphase <-GetGeneMultiple(gene = exemplary_genes,dataset = dataset_Brar,hdf5file = hdf5file_A_1,gffdf = gff_df,min_read_length = 10, asite_disp_length = asite_disp_length, scaling_factor = scaling_factor_A_1)

exemplary_genes_spore <-GetGeneMultiple(gene = exemplary_genes,dataset = dataset_Brar,hdf5file = hdf5file_SPORE_1,gffdf = gff_df,min_read_length = 10, asite_disp_length = asite_disp_length, scaling_factor = scaling_factor_SPORE_1)

# A tibble: 300 x 2
#     Pos Counts
#   <int>  <dbl>
# 1  -250      0
# 2  -249      0
# 3  -248      0
# 4  -247      0
# 5  -246      0
# 6  -245      0
# 7  -244      0
# 8  -243      0
# 9  -242      0
# 10  -241      0
# # … with 290 more rows

#GetGeneMeta
GetGeneMeta(gene = exemplary_genes,dataset = dataset_Brar,hdf5file = hdf5file_SPORE_1,gffdf = gff_df,min_read_length = 10, asite_disp_length = asite_disp_length, scaling_factor = scaling_factor_SPORE_1)

#join_GetGeneData
gene1<-GetGene(gene = exemplary_genes[1], dataset = dataset_Brar, hdf5file = hdf5file_SPORE_1,gffdf = gff_df,nnt_gene = nnt_gene,min_read_length = 10,asite_disp_length = asite_disp_length, snapdisp = 0L,scaling_factor = scaling_factor_SPORE_1)

gene2<- GetGene(gene = exemplary_genes[2], dataset = dataset_Brar, hdf5file = hdf5file_SPORE_1,gffdf = gff_df,nnt_gene = nnt_gene,min_read_length = 10,asite_disp_length = asite_disp_length,snapdisp = 0L,scaling_factor = scaling_factor_SPORE_1)

join_GetGeneData(gene1,gene2 )
# A tibble: 300 x 2
#     Pos Counts
#    <int>  <dbl>
# 1  -250      0
# 2  -249      0
# 3  -248      0
# 4  -247      0
# 5  -246      0
# 6  -245      0
# 7  -244      0
# 8  -243      0
# 9  -242      0
# 10  -241      0
# # … with 290 more rows

#CanVsNon
exempary_VEG_CanVsNon <-CanVsNon(exemplary_genes_vegetative)
exempary_A_CanVsNon <-CanVsNon(exemplary_genes_anaphase)
exempary_SPORE_CanVsNon <-CanVsNon(exemplary_genes_spore)
# A tibble: 6 x 3
# Gene    Count_TL Count_AUG
# <chr>      <dbl>     <dbl>
# 1 YCR012W  0           0    
# 2 YEL009C  0.0101      0.193
# 3 YGR094W  0           0    
# 4 YML065W  0.00506     0    
# 5 YOL130W  0.118       0    
# 6 YOR303W  0.223       0 

#TLvsAUGscatter
TLvsAUGscatter(exempary_genes_CanVsNon)

#TL_AUG_Regulation
TL_AUG_Regulation(exempary_SPORE_CanVsNon,exempary_A_CanVsNon,exempary_SPORE_CanVsNon )

#DistributionPlot
DistributionPlot(exemplary_genes_spore[[1]])

MultipleDistributionPlots(exemplary_genes_vegetative, exemplary_genes)

#TwoPlotsDistribution
TwoPlotsDistribution(exemplary_genes_vegetative[[1]], exemplary_genes_anaphase[[1]],gene ="YCR012W", cond1_name = "VEGETATIVE", "ANPHASE_1" )

#ThreePlotsDistribution
ThreePlotsDistribution(exemplary_genes_vegetative[[1]], exemplary_genes_anaphase[[1]], exemplary_genes_spore[[1]], cond1_name = "VEGETATIVE", "ANPHASE_1", "SPORE")

#sum_uAUG
sum_uAUG(datatibble = exemplary_genes_vegetative[[1]],uAUG_start = -250,uAUG_end = -5)
# # A tibble: 1 x 1
#     sum_read_counts_uAUG
#                    <dbl>
#   1               0.0132

#sum_uAUG_multiple
sum_uAUG_multiple(genes = exemplary_genes_vegetative,uAUG_start = -250, uAUG_end = -5)
# $YCR012W.sum_read_counts_uAUG
# [1] 0.01318392
# 
# $YEL009C.sum_read_counts_uAUG
# [1] 0.1379488
# 
# $YOR303W.sum_read_counts_uAUG
# [1] 0.02894031
# 
# $YOL130W.sum_read_counts_uAUG
# [1] 0.01575639
# 
# $YGR094W.sum_read_counts_uAUG
# [1] 0.0006431181
# 
# $YML065W.sum_read_counts_uAUG
# [1] 0.0003215591

#CDS_TL
exemplary_result <-CDS_TL(genes = exemplary_genes_vegetative,gene_names = exemplary_genes)
# # A tibble: 12 x 3
# genes     region    count
# <list>    <chr>     <dbl>
# 1 <chr [1]> TL     0.0132  
# 2 <chr [1]> TL     0.138   
# 3 <chr [1]> TL     0.0289  
# 4 <chr [1]> TL     0.0158  
# 5 <chr [1]> TL     0.000643
# 6 <chr [1]> TL     0.000322
# 7 <chr [1]> AUG    1.52    
# 8 <chr [1]> AUG    0.137   
# 9 <chr [1]> AUG    0       
# 10 <chr [1]> AUG    0       
# 11 <chr [1]> AUG    0.00525 
# 12 <chr [1]> AUG    0.0420 

#GetTLvsAUG
GetTLvsAUG(genes = exemplary_result,gene_names = exemplary_genes)

#UpstreamRegionsCount
exemplaru_regions_veg <-UpstreamRegionsCount(exemplary_genes_vegetative,gene_names = exemplary_genes)
# # A tibble: 36 x 3
# genes     region     read
# <list>    <chr>     <dbl>
# 1 <chr [1]> TL     0.0441  
# 2 <chr [1]> TL     0.142   
# 3 <chr [1]> TL     0.0284  
# 4 <chr [1]> TL     0.0154  
# 5 <chr [1]> TL     0.000630
# 6 <chr [1]> TL     0.00284 
# 7 <chr [1]> AUG    1.52    
# 8 <chr [1]> AUG    0.137   
# 9 <chr [1]> AUG    0       
# 10 <chr [1]> AUG   0       
# # … with 26 more rows

#Different_Regions_barplot
Different_Regions_barplot(data = exemplaru_regions_veg)

#plotted_each_barplot
plotted_each_barplot(exemplaru_regions_veg, exemplary_genes[[1]])

#all_genes_TL_AUG
all_genes_TL_AUG(exempary_SPORE_CanVsNon)

#TL_3_conditions_all
TL_3_conditions_all(condition1_TL = exempary_VEG_CanVsNon,condition2_TL =  exempary_A_CanVsNon, condition3_TL = exempary_SPORE_CanVsNon, name1 = "VEGETATIVE", "ANAPHASE I", "SPORE")

#compared_TL_single
compared_TL_single(condition1 = exempary_VEG_CanVsNon, condition2 = exempary_A_CanVsNon, condition3 = exempary_SPORE_CanVsNon, gene = "YCR012W")

#Efficiency_preparation
Efficiency_preparation(condition1 = exempary_VEG_CanVsNon,condition2 =  exempary_A_CanVsNon, condition3 = exempary_SPORE_CanVsNon)
# # A tibble: 18 x 4
# Gene    Condition Reads_TL Reads_AUG
# <chr>   <chr>        <dbl>     <dbl>
# 1 YCR012W VEG       0.0132     1.51   
# 2 YEL009C VEG       0.138      0.137  
# 3 YGR094W VEG       0.000322   0.00525
# 4 YML065W VEG       0.000322   0.0420 
# 5 YOL130W VEG       0.0158     0      
# 6 YOR303W VEG       0.0289     0      
# 7 YCR012W ANAPH_I   0          0      
# 8 YEL009C ANAPH_I   0.0222     0.214  
# 9 YGR094W ANAPH_I   0          0      
# 10 YML065W ANAPH_I  0.00262    0      
# 11 YOL130W ANAPH_I  0.0497     0      
# 12 YOR303W ANAPH_I  0.162      0      
# 13 YCR012W SPORE    0          0      
# 14 YEL009C SPORE    0.0101     0.193  
# 15 YGR094W SPORE    0          0      
# 16 YML065W SPORE    0.00506    0      
# 17 YOL130W SPORE    0.118      0      
# 18 YOR303W SPORE    0.223      0      

#Efficiency_AUG_TL
Efficiency_AUG_TL(condition1 = exempary_VEG_CanVsNon,condition2 =  exempary_A_CanVsNon, condition3 = exempary_SPORE_CanVsNon, gene = exemplary_genes[[2]], condition_name1 = "VEGETATIVE", "ANAPHASE I", "SPORE")
