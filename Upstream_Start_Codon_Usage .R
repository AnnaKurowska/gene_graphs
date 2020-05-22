
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
library(stringr)
library(ape)
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

# Takes value for start position of 5'LS for matrix creation (no negative values e.g.: 1)
GetTLstart <- function(gene, gffdf) {
  gffdf %>% 
    dplyr::filter(Name == gene, type=="LS5" ) %>% 
    dplyr::pull(start)
}

# Takes value for 3' end of the table (includes whole 5'LS +50 nt of CDS)
CDS3_end <- function(gene, gffdf) {
  gffdf %>% 
    dplyr::filter(Name == gene, type=="LS5") %>% 
    dplyr::pull(end) + nnt_gene
} 

# Creates a matrix with number of columns that start at 5' LS and finish at 50'th NT of  CDS
GetGeneDatamatrixTL <- function(gene, dataset, hdf5file, gffdf, nnt_gene) {
  data_mat_all <- GetGeneDatamatrix(gene, dataset, hdf5file)
  n_left5 <- GetTLstart(gene, gffdf) # column to start from (5'end)
  n_right3 <-CDS3_end(gene, gffdf)  # column to end with (3'end)
  data_mat_5start <- data_mat_all[, n_left5 : n_right3]
  # data_mat_5start <- tibble::as_tibble(data_mat_5start, .name_repair="minimal")
  return(data_mat_5start)
}

## A-site mapping function
CalcAsiteFixedOneLength <- function(reads_pos_length, min_read_length,
                                    read_length, asite_disp) {
  # Calculate read A-site using a fixed displacement for a single read length
  length_row_choose <- read_length - min_read_length + 1
  reads_pos_length[length_row_choose, ] %>%
    dplyr::lag(n = asite_disp, default = 0)
}


## A-site mapping function
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

## A-site mapping function
SnapToCodon <- function(x, left, right, snapdisp=0L) {
  # snap nucleotide-aligned reads to codon position
  #   x:     vector
  #   left:  integer for starting position, frame 0
  #   right: integer for ending position
  #   snapdisp: integer any additional displacement in the snapping
  RcppRoll::roll_suml(x[(left:right) + snapdisp], n=1L, by=1L, fill = NULL)
  
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
  # end = LS5_length(gene, gffdf)
  
  # will have to call x and gff_df the same!
  reads_pos_length <- GetGeneDatamatrixTL(gene, dataset,
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

##A-site mapping for multiple genes; output is an individual tibble for each gene
GetGeneMultiple <- function(gene,
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
                      scaling_factor = scaling_factor_SPORE_2)
  
  names(output) <- gene
  
  
  return(output)
}


## A-site mapping for multiple genes; output is one tibble with cumulative count of reads at each position
GetGeneMeta <- function(gene, dataset, hdf5file, gffdf,
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
                      scaling_factor = scaling_factor_SPORE_2)
  
  names(output) <- gene
  
  Counts <- output %>%
    Reduce("+", .) %>%
    dplyr::select(Counts)
  
  Pos <- output$YEL009C$Pos ### 
  
  Counts_Asite_mapped_all <-base::cbind(Pos, Counts) %>%
    tibble::as_tibble(.name_repair = "minimal")
  
  return(Counts_Asite_mapped_all)
  
}
## Combines 2 outputs from GetGeneCodonPosReads1dsnap or All_genesAmapped
join_reduced_data <-function(dataset1, dataset2){
  full_join(dataset1,dataset2, by = "Pos") %>%
    mutate(., Counts = Counts.x + Counts.y) %>%
    select(Pos, Counts)%>%
    return()
} 

#### functions which take output from A_mapped_genes as input

# Density in 5'LS (-250:-5) and AUG (-5: 10)
CanVsNon <- function(dataset) {
  # Genes_names <-c(gene)
  
  LS <- bind_rows(dataset, .id = "Gene") %>%
    dplyr::filter( Pos >= -250, Pos < -5) %>%
    dplyr::group_by(Gene) %>%
    dplyr::summarise(Count_TL = sum(Counts)/(length(-250:-5)-1))
  
  AUG <- bind_rows(dataset, .id = "Gene") %>%
    dplyr::filter( Pos >= -5, Pos < 10) %>%
    dplyr::group_by(Gene) %>%
    dplyr::summarise(Count_AUG = sum(Counts)/(length(-5:10)-1))
  
  LS_AUG <-dplyr::full_join(LS, AUG, by = "Gene") 
  
}

###a scatter plot showing the LS/AUG ratio, input: CanVsNon output
scatterAUG_LS_ratio <-function(CanNonData){
  ggplot(CanNonData, aes(x = Counts_AUG, y = Counts_LS5)) +
    geom_point() +
    # scale_y_continuous(expand = c(0,0)) +
    scale_y_log10() +
    scale_x_log10() +
    geom_text_repel(aes(label = Genes)) +
    labs(x= "AUG (CPMB)", y = "TL (CPMB)") +
    theme_classic()
  
}

##### Same gene across different stages, input: output of CanVsNon data for each stage
LS5_AUG_Regulation <-function(CanNonData1, CanNonData2, CanNonData3){
  
  comparing_LS_test <- full_join(CanNonData1, CanNonData2, by = "Gene") %>%
    full_join(CanNonData3, by = "Gene") %>%
    set_colnames(c("Gene", "LS_WT_NONE", "AUG_WT_NONE", "LS_WT_CHX","AUG_WT_CHX", "LS_WT_3AT", "AUG_WT_3AT" )) 
  
  comparing_LS_test <-comparing_LS_test%>%
    select("Gene","LS_WT_NONE","LS_WT_CHX","LS_WT_3AT") %>%
    gather(key = "Condition", value = "Reads", -Gene)
  
  comparing_LS_test$Condition <- factor(comparing_LS_test$Condition, levels=unique(comparing_LS_test$Condition))  
  
  ggplot(comparing_LS_test, aes(x = Condition, group = Gene, y= Reads, color = Gene)) +
    geom_point() +
    geom_line() +
    scale_y_continuous(expand = c(0,0)) +
    geom_text_repel(aes(label = round(Reads, digits = 3)))
  
}
## distribution plot, input: output from GetGeneCodonPosReads1dsnap or All_genesAmapped

plotting_TL<- function(input_data) {
  
  plotted_LS <- ggplot2::ggplot(input_data) +
    geom_density(aes(x=Pos, y=Counts), stat="identity") +
    scale_y_continuous(expand = c(0,0)) +
    scale_x_continuous(breaks = c(-250,-200,-150,-100,-60,-30,0,30,50), labels=c("-250","200","-150","100", "-60","-30","0\nAUG",30,50)) +
    labs( x = "Position", y = "CPM") +
    theme_classic() 
  return(plotted_LS)
}



## Multiple distribution plots, input: output from A_mapped_genes, GetGeneCodonPosReads1dsnap or All_genesAmapped

plotting_multiple <- function(input_data, gene_names) {
  purrr :: map(input_data, plotting_TL) %>%
    ggarrange(plotlist = ., labels = gene_names) %>%
    return()
}

## distribution plot of 2 samples, input: GetGeneCodonPosReads1dsnap,A_mapped_genes (1 gene only) or All_genesAmapped

two_in_one_plots <- function(sample1, sample2, gene, cond1_name, cond2_name) {
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

## distribution plot of 3 samples, input: GetGeneCodonPosReads1dsnap,A_mapped_genes (1 gene only) or All_genesAmapped

three_in_one_plots <- function(sample1, sample2, sample3, cond1_name, cond2_name, cond3_name, gene) {
  ggplot2::ggplot(sample1,aes(x=Pos, y=Counts, alpha = 0.5)) +
    geom_density(data = sample1, stat="identity", aes(color = cond1_name, fill = cond1_name )) +
    geom_density(data = sample2, stat="identity", aes(color = cond2_name, fill = cond2_name)) +
    geom_density(data = sample3, stat="identity", aes(color = cond3_name, fill =cond3_name )) +
    scale_x_continuous(breaks = c(-250,-200,-150,-100,-60,-30,0,30,50), labels=c("-250","200","-150","100", "-60","-30","0\nAUG",30,50)) +
    scale_y_continuous(expand = c(0,0)) +
    labs(y= "CPM", x = "Position") +
    theme_classic() +
    scale_color_manual(values = c(VEG = "red", ANAPH_I = "blue", SPORE = "black" ))
}


## sum_uAUG and sum_uAUG_mulitple used in CDS_TL and all_regions_together

##Sums read counts and normalizes them for density at 5'LS and AUG 
sum_uAUG <- function(datatibble, uAUG_start, uAUG_end){
  datatibble %>%
    dplyr::filter(Pos >= uAUG_start, Pos <= uAUG_end) %>%
    dplyr::summarise(sum_read_counts_uAUG=sum(Counts)/(length(uAUG_start:uAUG_end)-1)) 
}

##Iteration of sum_uAUG over many genes
sum_uAUG_multiple<- function(genes, uAUG_start, uAUG_end) {
  sapply(genes, sum_uAUG, uAUG_start, uAUG_end) %>%
    return()
}

## 

CDS_TL <- function(genes, gene_names) {
  table_genes <-as.character((paste0(gene_names)))  #why does it make it a factor?
  
  
  region1 <-sum_uAUG_multiple(genes, uAUG_start = -250, uAUG_end = -5)
  region2 <- sum_uAUG_multiple(genes, uAUG_start = -5, uAUG_end = 10)
  
  
  result <-base::cbind(table_genes, region1, region2) %>%
    as_tibble() %>%
    set_colnames(c("genes","5'TL", "AUG"))
  
  result_new <-result %>%
    tidyr::gather( key = "region", value = "count", -genes)
  result_new$count <- as.numeric(result_new$count)
  return(result_new)
}

## bar plot showing AUG and TL for individual genes, requires separate output from CDS_TL
efficiency_barplot_each_gene <- function(CDS_TL_AUG_data, gene){
  positions <- c( "5'LS", "AUG")
  value <-filter(data, genes == gene)
  # 
  ggplot(CDS_TL_AUG_data, aes(x = region, y = count, fill = region)) +
    geom_col(show.legend = FALSE) +
    theme_classic() +
    labs(y= "CPMB", x = "Region (mRNA)") +
    # ggtitle(paste0(gene)) +
    scale_y_continuous(expand = c(0,0)) +
    scale_x_discrete(limits = positions)
}

# efficiency bar plot (TL and AUG) for all input genes, includes the CDS_TL step; input: GetGeneCodonPosReads1dsnap,A_mapped_genes (1 gene only) or All_genesAmapped
efficiency_barplot <- function(genes, gene_names) {
  CDS_TL_data <-CDS_TL(genes, gene_names) %>%
    tidyr::gather( key = "region", value = "count", -genes)
  CDS_TL_data$count <- as.numeric(CDS_TL_data$count)
  
  positions <- c("5'LS", "AUG")
  
  ggplot(CDS_TL_data, aes(x=region, y = count, fill = region )) +
    geom_col(show.legend = FALSE) +
    theme_classic() +
    labs(y= "CPMB", x = "Region (mRNA)") +
    scale_y_continuous(expand = c(0,0)) +
    scale_x_discrete(limits = positions) 
  
}


#Creates a table with read density at different regions, input:GetGeneCodonPosReads1dsnap,A_mapped_genes (1 gene only) or All_genesAmapped

all_regions_together <- function(genes, gene_names) {
  table_genes <-as.character((paste0(gene_names)))  #why does it make it a factor?
  
  region1 <-sum_uAUG_multiple(genes, uAUG_start = -250, uAUG_end = 0)
  region2 <- sum_uAUG_multiple(genes, uAUG_start = -5, uAUG_end = 10)
  region3 <-sum_uAUG_multiple(genes, uAUG_start = -15, uAUG_end = 0)
  region4 <-sum_uAUG_multiple(genes, uAUG_start = -60, uAUG_end = 0)
  region5 <-sum_uAUG_multiple(genes, uAUG_start = -120, uAUG_end = -60)
  region6  <-sum_uAUG_multiple(genes, uAUG_start = -250, uAUG_end = -120)
  
  
  created_table <-base::cbind(table_genes, region1, region2, region3, region4, region5, region6) %>%
    as_tibble() %>%
    magrittr::set_colnames(c("genes","TL", "AUG", "-15:-0", "-60:0", "-120:-60", "-250:-60")) %>%
    tidyr::gather( key = "region", value = "read", -genes)
  created_table$read <- as.numeric(created_table$read)
  
  return(created_table)
}


##Bar plots showing read density at each region for all selected genes (all or multiple), input: output of all_regions_together
plotted_barplots <- function(data) {
  
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


## Comparison of 5'LS to AUG for all or multiple selected genes in 1 condition, input: Output of CanVsNon
all_genes_LS_AUG <-function(condition1_LS){
  LS <- sum(condition1_LS$Counts_LS5)
  AUG <- sum(condition1_LS$Counts_AUG)
  
  Condition <- c("TL", "AUG")
  for_plot_LS <- rbind(LS, AUG) %>%
    as_tibble %>%
    set_colnames("LS5") %>%
    cbind(Condition)
  
  plot_LS <-ggplot(for_plot_LS) +
    geom_bar(aes(x=Condition, y= LS5, fill = Condition), stat="identity", show.legend = FALSE ) +
    # scale_x_discrete(labels=c("WT_3AT","WT_CHX", "WT_NONE")) +
    scale_y_continuous(expand = c(0,0)) +
    theme_classic() +
    labs(x= "Region", y = "TL (CPMB) ")
  return(plot_LS)
}


### Comparison of 5'LS for all genes between 3 conditions, input: output of CanVsNon 

LS5_3_conditions_all <-function(condition1_LS, condition2_LS, condition3_LS, name1, name2, name3){
  condition1 <- sum(condition1_LS$Counts_LS5)
  condition2 <- sum(condition2_LS$Counts_LS5)
  condition3 <- sum(condition3_LS$Counts_LS5)
  
  Condition <- paste0(c(name1, name2, name3))
  # Condition <- c("WT NONE", "WT CHX", "WT 3AT")
  for_plot_LS <- rbind(condition1, condition2, condition3) %>%
    as_tibble %>%
    set_colnames("LS5") %>%
    cbind(Condition)
  
  for_plot_LS$Condition <- factor(for_plot_LS$Condition, levels=unique(for_plot_LS$Condition))  
  
  plot_LS <-ggplot(for_plot_LS) +
    geom_bar(aes(x=Condition, y= LS5, fill = Condition), stat="identity", show.legend = FALSE ) +
    # scale_x_discrete(labels=c("WT_3AT","WT_CHX", "WT_NONE")) +
    scale_y_continuous(expand = c(0,0)) +
    theme_classic() +
    labs(x= "Stage", y = "TL (CPMB)")
  return(plot_LS)
}


### Barplot Comparison of 5'LS for 1 gene between 3 conditions, input: output of CanVsNon  

compared_LS_single<-function(condition1, condition2, condition3, gene){
  
  search_for_TL <- function(condition, gene){
    condition %>%
      filter(Gene == gene ) %>%
      select(Counts_LS5) %>%
      pull
  }
  
  LS1 <-search_for_TL(condition1, gene)
  
  LS2 <-search_for_TL(condition2, gene)
  
  LS3 <-search_for_TL(condition3, gene)
  
  Condition <- c("WT NONE", "WT CHX", "WT 3AT")
  
  for_plot <- rbind(LS1,LS2,LS3) %>%
    set_colnames("LS5") %>%
    cbind(Condition) %>%
    as_tibble()
  for_plot$LS5 <-as.numeric(for_plot$LS5)
  
  plot_LS <-ggplot(for_plot) +
    geom_bar(aes(x=Condition, y= LS5, fill = Condition), stat="identity",show.legend = FALSE ) +
    # scale_x_discrete(labels=c("META_I_1", "META_II_1", "PRE_ENTRY_1")) +
    scale_y_continuous(expand = c(0,0)) +
    theme_classic() +
    labs(x= "Condition", y = "CPMB 5'LS") +
    ggtitle(paste0(gene))
  
  return(plot_LS)
  
}


##function used within LS_AUG_plot_correlation and LS_AUG_correlation; input: 

finding_correlation_preparation<- function(gene_VEG, gene_Anaph, gene_Spore){
  
  comparing_genes_stages <- full_join(gene_VEG, gene_Anaph, by = "Gene") %>%
    full_join(gene_Spore, by = "Gene") %>%
    set_colnames(c("Gene", "VEG_LS", "VEG_AUG", "ANAPH_I_LS","ANAPH_I_AUG", "SPORE_LS",
                   "SPORE_AUG" )) 
  
  comparing_genes_stages_TL <-comparing_genes_stages %>%
    select("Gene","VEG_LS","ANAPH_I_LS","SPORE_LS") %>%
    set_colnames(c("Gene", "VEG", "ANAPH_I", "SPORE")) %>%
    gather(key = "Condition", value = "Reads_LS", -Gene) 
  
  comparing_genes_stages_AUG <-comparing_genes_stages %>%
    select("Gene","VEG_AUG","ANAPH_I_AUG","SPORE_AUG") %>%
    set_colnames(c("Gene", "VEG", "ANAPH_I", "SPORE")) %>%
    gather(key = "Condition", value = "Reads_AUG", -Gene)
  
  
  combined_LS_AUG <- full_join(comparing_genes_stages_TL, 
                                comparing_genes_stages_AUG, by = "Condition") %>%
    select(Gene.x, Condition, Reads_LS, Reads_AUG)
  
  return(combined_LS_AUG)
}


### outputs a figure with 5'LS and AUG densities across 3 conditions, input: output of CanVsNon
LS_AUG_plot_correlation <- function(gene_VEG, gene_Anaph, gene_Spore, condition_name1,condition_name2, condition_name3 ){
  
  prepared_table <- finding_correlation_preparation(gene_VEG, gene_Anaph, gene_Spore)
  
  prepared_table$Condition <- factor(prepared_table$Condition, levels=unique(prepared_table$Condition))
  
  plot <-ggplot(prepared_table) +
    geom_point(aes(x = Condition, group = Gene.x, y= Reads_AUG, color = "AUG" )) +
    geom_line(aes(x = Condition, group = Gene.x, y= Reads_AUG, color = "AUG")) +
    geom_point(aes(x = Condition, group = Gene.x, y= Reads_LS, color = "TL")) +
    geom_line(aes(x = Condition, group = Gene.x, y= Reads_LS, color = "TL")) +
    scale_y_log10() +
    labs(y = "CPMB", x = "Condition") +
    guides(color=guide_legend(NULL)) +
    scale_x_discrete(labels= c(paste0(condition_name1), paste0(condition_name2), paste0(condition_name3))) +
    theme_classic()
  return(plot)
  
}
##calculates correlation between the 5'LS and AUG densities across 3 conditions, input: output of CanVsNon
LS_AUG_correlation <- function(gene_VEG, gene_Anaph, gene_Spore){
  
  prepared_table <- finding_correlation_preparation(gene_VEG, gene_Anaph, gene_Spore)
  
  cor.test(prepared_table$Reads_AUG, prepared_table$Reads_LS, method=c("pearson")) %>%
    return()
}



# ################### A SITE MAPPING #########

# 
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
  n_right3 <- LS5_length(gene, gffdf) + nnt_gene # column to end with (3'end)
  data_mat_5start <- data_mat_all[, n_left5 : n_right3]
  data_mat_5start <- tibble::as_tibble(data_mat_5start, .name_repair="minimal")
  return(data_mat_5start)
}

LS5_table <- function(gene) {
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


## iteration of LS5_table over multiple genes
final_function_table <- function(genes){
  table <-purrr::map(genes, LS5_table)
  names(table) <- genes
  return(table)
  
}

plotting_meta_analysis<- function(input_data) {
  
  plotted_LS <- ggplot(input_data) +
    geom_density(aes(x=Pos, y=Counts), stat="identity") +
    scale_x_continuous(limits = c(-250,50), expand = c(0, 0)) +
    scale_y_continuous(expand = c(0,0)) +
    labs(y= "Read count", x = "Position") +
    ggtitle(paste0("Ribosome footprint density of all genes: WT 3-AT")) +
    theme_classic() %>%
    return()
}

################################### Set up ######################################


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

############################ Brar et al. (2012)

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
#Initial set ups
dataset_G2014 <- "G-Sc_2014"
dataset_Brar <- "B-Sc_2012"
genes <- c("YPR036W-A", "YEL009C", "YOR303W", "YGR094W", "YOL130W", "YOR061W", "YAL040C", "YKL109W", "YIL144W", "YGR037C", "YJL106W", "YBR160W", "YDR172W", "YML065W", "YBR257W", "YOL104C", "YJR094C","YOL125W")
gene_names <- rhdf5::hTL(hdf5file_none, recursive = 1)$name

nnt_gene<- 50
startpos <-250
startlen <- 10
temporary_length <- 20 #for zooming in bit
min_read_count <- 10
orf_gff_file <- "G-Sc_2014/input/yeast_CDS_w_250utrs.gff3"

gff_df <- readGFFAsDf(orf_gff_file) 

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
