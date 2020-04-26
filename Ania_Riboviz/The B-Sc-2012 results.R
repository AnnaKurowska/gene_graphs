####################################### SAMPLE PREPARATION###################################

library(here)
here()
# META_I_1
META_I_1 <- "B-Sc_2012/input/output/META_I_1/META_I_1.h5"
hd_file_M_I_1 <- META_I_1
hdf5file_M_I_1 <- rhdf5::H5Fopen(hd_file_M_I_1)

#META_I_2

META_I_2 <- "B-Sc_2012/input/output/META_I_2/META_I_2.h5"
hd_file_M_I_2 <- META_I_2
hdf5file_M_I_2 <- rhdf5::H5Fopen(hd_file_M_I_2)

#META_II_1

META_II_1 <- "B-Sc_2012/input/output/META_II_1/META_II_1.h5"
hd_file_M_II_1 <- META_II_1
hdf5file_M_II_1 <- rhdf5::H5Fopen(hd_file_M_II_1)

#META_II_2

META_II_2 <- "B-Sc_2012/input/output/META_II_2/META_II_2.h5"
hd_file_M_II_2 <- META_II_2
hdf5file_M_II_2 <- rhdf5::H5Fopen(hd_file_M_II_2)

#PRE_ENTRY1

PRE_ENTRY_1 <- "B-Sc_2012/input/output/PRE_ENTRY_1/PRE_ENTRY_1.h5"
hd_file_PRE_ENTRY_1 <- PRE_ENTRY_1
hdf5file_PRE_ENTRY_1 <- rhdf5::H5Fopen(hd_file_PRE_ENTRY_1)

#PRE_ENTRY2

PRE_ENTRY_2 <- "B-Sc_2012/input/output/PRE_ENTRY_2/PRE_ENTRY_2.h5"
hd_file_PRE_ENTRY_2 <- PRE_ENTRY_2
hdf5file_PRE_ENTRY_2 <- rhdf5::H5Fopen(hd_file_PRE_ENTRY_2)


dataset_Brar <- "B-Sc_2012"
test_orfs <- c("YCR012W","YEL009C","YOR303W","YOL130W","YGR094W", "YML065W")
gene_names <- rhdf5::h5ls(hdf5file_M_I_1, recursive = 1)$name
gff_df <- readGFFAsDf(orf_gff_file) 

 
############################################A SITE MAPPING##################################
# meta 1
genes <- c("YPR036W-A", "YEL009C", "YOR303W", "YGR094W", "YOL130W", "YOR061W", "YAL040C", "YKL109W", "YIL144W", "YGR037C", "YJL106W", "YBR160W", "YDR172W", "YML065W", "YBR257W", "YOL104C", "YJR094C")

#sevR = several genes reduced
#sev = several genes (not reduced)
#allR = all genes reduced

######1 META_I_1 1) several genes together 2) several genes separately 3) all genes together ######

META_I_1_sevR <-All_genesAmapped(gene = genes, dataset = dataset_Brar, hdf5file = hdf5file_M_I_1, gffdf = gff_df, min_read_length = 10)
#if it gives an error, that's because of the output$YEL009C$Pos

META_I_1_sev <- A_mapped_genes(genes, dataset = dataset_Brar, hdf5file = hdf5file_M_I_1, gffdf = gff_df, min_read_length = 10)

META_I_1_sev_5UTR <- CanVsNon(META_I_1_sev)

META_I_1_allR <-All_genesAmapped(gene = gene_names, dataset = dataset_Brar, hdf5file = hdf5file_M_I_1, gffdf = gff_df, min_read_length = 10)

META_I_1_allR_plot <-plotting_5UTR(META_I_1_allR)
ggsave("META_I_1_allR_plot", device = "jpg")

META_I_1_sev_5UTR_scatter <- ggplot(META_I_1_sev_5UTR, aes(x = Counts_AUG, y = Counts_UTR5)) +
  geom_point() + 
  geom_text(aes(label = Gene))
ggsave("META_I_1_sev_5UTR_scatter", device = "jpg")
# scale_y_continuous(expand = c(0,0))

META_I_1_sev_regions <- all_regions_together(META_I_1_sev,gene_names = genes) 
plotted_barplots(META_I_1_sev_regions) 
#Error in is.finite(x) : default method not implemented for type 'list'

plotting_multiple(input_data = META_I_1_sev, genes)
ggsave("META_I_1_all_genes", device = "jpg")

META_I_1_sev_UTR_AUG_plot <-CDS_5UTR(genes = META_I_1_sev, gene_names = genes) %>%
  tidyr::gather( key = "region", value = "count", -genes)
META_I_1_sev_UTR_AUG_plot$count <- as.numeric(META_I_1_sev_UTR_AUG_plot$count) 
META_I_1_sev_UTR_AUG_plot <- efficiency_barplot(META_I_1_sev_UTR_AUG_plot)
ggsave("META_I_1_sev_UTR_AUG_plot", device = "jpg")

efficiency_barplot_each_gene(together, "YJR094C") ##!!! ughhh

########2 META_I_2 1) several genes together 2) several genes separately 3) all genes together##### 
META_I_2_sevR <-All_genesAmapped(gene = genes, dataset = dataset_Brar, hdf5file = hdf5file_M_I_2, gffdf = gff_df, min_read_length = 10)
#if it gives an error, that's because of the output$YEL009C$Pos

META_I_2_sev <- A_mapped_genes(genes, dataset = dataset_Brar, hdf5file = hdf5file_M_I_2, gffdf = gff_df, min_read_length = 10)

META_I_2_allR <-All_genesAmapped(gene = gene_names, dataset = dataset_Brar, hdf5file = hdf5file_M_I_2, gffdf = gff_df, min_read_length = 10)

META_I_2_allR_plot <-plotting_5UTR(META_I_2_allR)
ggsave("META_I_2_allR_plot", device = "jpg")

META_I_2_sev_5UTR <- CanVsNon(META_I_2_sev)

META_I_2_sev_5UTR_scatter <- ggplot(META_I_2_sev_5UTR, aes(x = Counts_AUG, y = Counts_UTR5)) +
  geom_point() +
  geom_text(aes(label = Gene))
ggsave("META_I_2_sev_5UTR_scatter", device = "jpg")

# scale_y_continuous(expand = c(0,0))

META_I_2_sev_regions <- all_regions_together(META_I_2_sev,gene_names = genes)

plotting_multiple(input_data = META_I_2_sev, genes)
ggsave("META_I_2_all_genes", device = "jpg")

META_I_2_sev_UTR_AUG_plot <-CDS_5UTR(genes = META_I_2_sev, gene_names = genes) %>%
  tidyr::gather( key = "region", value = "count", -genes)
META_I_2_sev_UTR_AUG_plot$count <- as.numeric(META_I_2_sev_UTR_AUG_plot$count) 
META_I_2_sev_UTR_AUG_plot <- efficiency_barplot(META_I_2_sev_UTR_AUG_plot)
ggsave("META_I_2_sev_UTR_AUG_plot", device = "jpg")

efficiency_barplot_each_gene(together, "YJR094C")

########3 META_II_1 1) several genes together 2) several genes separately 3) all genes together#### 
META_II_1_sevR <-All_genesAmapped(gene = genes, dataset = dataset_Brar, hdf5file = hdf5file_M_II_1, gffdf = gff_df, min_read_length = 10)
#if it gives an error, that's because of the output$YEL009C$Pos

META_II_1_sev <- A_mapped_genes(genes, dataset = dataset_Brar, hdf5file = hdf5file_M_II_1, gffdf = gff_df, min_read_length = 10)

META_II_1_allR <-All_genesAmapped(gene = gene_names, dataset = dataset_Brar, hdf5file = hdf5file_M_II_1, gffdf = gff_df, min_read_length = 10)

META_II_1_allR_plot <-plotting_5UTR(META_II_1_allR)
ggsave("META_II_1_allR_plot", device = "jpg")

META_II_1_sev_5UTR <- CanVsNon(META_II_1_sev)

META_II_1_sev_5UTR_scatter <- ggplot(META_II_1_sev_5UTR, aes(x = Counts_AUG, y = Counts_UTR5)) +
  geom_point() +
  geom_text(aes(label = Gene))
ggsave("META_II_1_sev_5UTR_scatter", device = "jpg")
# scale_y_continuous(expand = c(0,0))

META_II_1_sev_regions <- all_regions_together(META_II_1_sev,gene_names = genes)

plotting_multiple(input_data = META_II_1_sev, genes)
ggsave("META_II_1_all_genes", device = "jpg")

META_II_1_sev_UTR_AUG_plot <-CDS_5UTR(genes = META_II_1_sev, gene_names = genes) %>%
  tidyr::gather( key = "region", value = "count", -genes)
META_II_1_sev_UTR_AUG_plot$count <- as.numeric(META_II_1_sev_UTR_AUG_plot$count) 
META_II_1_sev_UTR_AUG_plot <- efficiency_barplot(META_II_1_sev_UTR_AUG_plot)
ggsave("META_II_1_sev_UTR_AUG_plot", device = "jpg")

efficiency_barplot_each_gene(together, "YJR094C")

####4 META_II_2 1) several genes together 2) several genes separately 3) all genes together########

META_II_2_sevR <-All_genesAmapped(gene = genes, dataset = dataset_Brar, hdf5file = hdf5file_M_II_2, gffdf = gff_df, min_read_length = 10)
#if it gives an error, that's because of the output$YEL009C$Pos

META_II_2_sev <- A_mapped_genes(genes, dataset = dataset_Brar, hdf5file = hdf5file_M_II_2, gffdf = gff_df, min_read_length = 10)

META_II_2_allR <-All_genesAmapped(gene = gene_names, dataset = dataset_Brar, hdf5file = hdf5file_M_II_2, gffdf = gff_df, min_read_length = 10)

META_II_2_allR_plot <-plotting_5UTR(META_II_2_allR)
ggsave("META_II_2_allR_plot", device = "jpg")

META_II_2_sev_5UTR <- CanVsNon(META_II_2_sev)

META_II_2_sev_5UTR_scatter <- ggplot(META_II_2_sev_5UTR, aes(x = Counts_AUG, y = Counts_UTR5)) +
  geom_point() +
  geom_text(aes(label = Gene))
ggsave("META_II_2_sev_5UTR_scatter", device = "jpg")
# scale_y_continuous(expand = c(0,0))

META_II_2_sev_regions <- all_regions_together(META_II_2_sev,gene_names = genes)

plotting_multiple(input_data = META_II_2_sev, genes)
ggsave("META_II_2_all_genes", device = "jpg")

META_II_2_sev_UTR_AUG_plot <-CDS_5UTR(genes = META_II_2_sev, gene_names = genes) %>%
  tidyr::gather( key = "region", value = "count", -genes)
META_II_2_sev_UTR_AUG_plot$count <- as.numeric(META_II_2_sev_UTR_AUG_plot$count) 
META_II_2_sev_UTR_AUG_plot <- efficiency_barplot(META_II_2_sev_UTR_AUG_plot)
ggsave("META_II_2_sev_UTR_AUG_plot", device = "jpg")

efficiency_barplot_each_gene(together, "YJR094C")

######5 PRE_ENTRY_1 1) several genes together 2) several genes separately 3) all genes together####

PRE_ENTRY_1_sevR <-All_genesAmapped(gene = genes, dataset = dataset_Brar, hdf5file = hdf5file_PRE_ENTRY_1, gffdf = gff_df, min_read_length = 10)
#if it gives an error, that's because of the output$YEL009C$Pos

PRE_ENTRY_1_sev <- A_mapped_genes(genes, dataset = dataset_Brar, hdf5file = hdf5file_PRE_ENTRY_1, gffdf = gff_df, min_read_length = 10)

PRE_ENTRY_1_allR <-All_genesAmapped(gene = gene_names, dataset = dataset_Brar, hdf5file = hdf5file_PRE_ENTRY_1, gffdf = gff_df, min_read_length = 10)

PRE_ENTRY_1_allR_plot <-plotting_5UTR(PRE_ENTRY_1_allR)
ggsave("PRE_ENTRY_1_allR_plot", device = "jpg")

PRE_ENTRY_1_sev_5UTR <- CanVsNon(PRE_ENTRY_1_sev)

PRE_ENTRY_1_sev_5UTR_scatter <- ggplot(PRE_ENTRY_1_sev_5UTR, aes(x = Counts_AUG, y = Counts_UTR5)) +
  geom_point() +
  geom_text(aes(label = Gene))
# scale_y_continuous(expand = c(0,0))
ggsave("PRE_ENTRY_1_sev_5UTR_scatter", device = "jpg")

PRE_ENTRY_1_sev_regions <- all_regions_together(PRE_ENTRY_1_sev,gene_names = genes)

plotting_multiple(input_data = PRE_ENTRY_1_sev, genes)
ggsave("PRE_ENTRY_1_all_genes", device = "jpg")

PRE_ENTRY_1_sev_UTR_AUG_plot <-CDS_5UTR(genes = PRE_ENTRY_1_sev, gene_names = genes) %>%
  tidyr::gather( key = "region", value = "count", -genes)
PRE_ENTRY_1_sev_UTR_AUG_plot$count <- as.numeric(PRE_ENTRY_1_sev_UTR_AUG_plot$count) 
PRE_ENTRY_1_sev_UTR_AUG_plot <- efficiency_barplot(PRE_ENTRY_1_sev_UTR_AUG_plot)
ggsave("PRE_ENTRY_1_sev_UTR_AUG_plot", device = "jpg")

efficiency_barplot_each_gene(together, "YJR094C")

####6 PRE_ENTRY_2 1) several genes together 2) several genes separately 3) all genes together#####

PRE_ENTRY_2_sevR <-All_genesAmapped(gene = genes, dataset = dataset_Brar, hdf5file = hdf5file_PRE_ENTRY_2, gffdf = gff_df, min_read_length = 10)
#if it gives an error, that's because of the output$YEL009C$Pos

PRE_ENTRY_2_sev <- A_mapped_genes(genes, dataset = dataset_Brar, hdf5file = hdf5file_PRE_ENTRY_2, gffdf = gff_df, min_read_length = 10)

PRE_ENTRY_2_allR <-All_genesAmapped(gene = gene_names, dataset = dataset_Brar, hdf5file = hdf5file_PRE_ENTRY_2, gffdf = gff_df, min_read_length = 10)

PRE_ENTRY_2_allR_plot <-plotting_5UTR(PRE_ENTRY_2_allR)
ggsave("PRE_ENTRY_2_allR_plot", device = "jpg")

PRE_ENTRY_2_sev_5UTR <- CanVsNon(PRE_ENTRY_2_sev)

PRE_ENTRY_2_sev_5UTR_scatter <-ggplot(PRE_ENTRY_2_sev_5UTR, aes(x = Counts_AUG, y = Counts_UTR5)) +
  geom_point() +
  geom_text(aes(label = Gene))
# scale_y_continuous(expand = c(0,0))
ggsave("PRE_ENTRY_2_sev_5UTR_scatter", device = "jpg")

PRE_ENTRY_2_sev_regions <- all_regions_together(PRE_ENTRY_2_sev,gene_names = genes)

PRE_ENTRY_2_sev_plot <-plotting_multiple(input_data = PRE_ENTRY_2_sev, genes)
ggsave("PRE_ENTRY_2_all_genes", device = "jpg")

PRE_ENTRY_2_sev_UTR_AUG_plot <-CDS_5UTR(genes = PRE_ENTRY_2_sev, gene_names = genes) %>%
  tidyr::gather( key = "region", value = "count", -genes)
PRE_ENTRY_2_sev_UTR_AUG_plot$count <- as.numeric(PRE_ENTRY_2_sev_UTR_AUG_plot$count) 
PRE_ENTRY_2_sev_UTR_AUG_plot <- efficiency_barplot(PRE_ENTRY_2_sev_UTR_AUG_plot)
ggsave("PRE_ENTRY_2_sev_UTR_AUG_plot", device = "jpg")

efficiency_barplot_each_gene(together, "YJR094C")

####################################### COMPARISONS ###############################
### 5UTR comparison for all genes (with META_I_1_sev_5UTR result)
compared_UTR <- UTR5_3_conditions_all(PRE_ENTRY_1_sev_5UTR, META_I_1_sev_5UTR, META_II_1_sev_5UTR)
ggsave("compared_UTR", device = "jpg")
 #sprawdzone ze kolejosc jest w porzadku 

### 5UTR comparison for one gene genes (with META_I_1_sev_5UTR result)
compared_UTR_YAL040C <- compared_UTR_single(META_I_1_sev_5UTR, META_I_2_sev_5UTR, PRE_ENTRY_1_sev_5UTR, "YML065W") 
ggsave("compared_UTR_YAL040C", device = "jpg")

### all genes plotted between 3 conditions

plot_BRAR <- ggarrange(META_I_1_allR_plot, META_I_2_allR_plot, META_II_1_allR_plot, META_II_2_allR_plot, PRE_ENTRY_1_allR_plot, PRE_ENTRY_2_allR_plot, labels = c("PRE-MEIOTIC ENTRY BioRep 1 (0h)","PRE-MEIOTIC ENTRY BioRep 2 (0h)", "METAPHASE I BioRep 1 (xh)", "METAPHASE I BioRep 2 (xh)", "METAPHASE II BioRep 1 (xh)", "METAPHASE II BioRep (xh)"), ncol = 1, nrow = 6 )
ggsave("plot_BRAR_all_conditions", device = "jpg")

#same gene across several conditions
#let's look at gene ORC1 YML065W

something <-GetGeneCodonPosReads1dsnap(gene = "YGR037C", dataset = dataset_Brar, hdf5file = hdf5file_M_II_1, gffdf = gff_df, nnt_gene = nnt_gene, min_read_length = 10, asite_disp_length = asite_disp_length)

something_plot <-plotting_5UTR(something)


################ quantitative comparisons in 5UTRs between conditions############

PRE_ENTRY_1_sev_5UTR
# A tibble: 17 x 3
#     Gene      Counts_UTR5 Counts_AUG
#     <chr>           <dbl>      <dbl>
# 1 YAL040C             0          7
# 2 YBR160W             1          0
# 3 YBR257W             0          2
# 4 YDR172W             0          4
# 5 YEL009C            30          0
# 6 YGR037C             5         70

comparing_UTR <- full_join(PRE_ENTRY_1_sev_5UTR, META_I_1_sev_5UTR, by = "Gene") %>%
  full_join(META_II_1_sev_5UTR, by = "Gene") %>%
  set_colnames(c("Gene", "UTR_PRE_ENTRY_1", "AUG_PRE_ENTRY_1", "UTR_META_I_1","AUG_META_I_1", "UTR_META_II_1", "AUG_META_II_1" )) 

why <-comparing_UTR %>%
  select("Gene","UTR_PRE_ENTRY_1","UTR_META_I_1","UTR_META_II_1") %>%
  gather(key = "Condition", value = "Reads", -Gene) %>%
  group_by(Gene)

ggplot(why, aes(x = Condition, y= Reads, color = Gene)) +
  geom_point() 
##no clue how to make these points connected in order to get one line & how to make each UTR column substract itself automatically from the other value


PRE_ENTRY_2_sev_regions
#
#    genes   region count
# 1	YPR036W-A	5'UTR	  17
# 2	YEL009C	  5'UTR	  30
# 3	YOR303W	  5'UTR	  209
# 4	YGR094W	  5'UTR	   0
# 5	YOL130W	  5'UTR 	 6
# 6	YOR061W	  5'UTR	   2
# 7	YAL040C	  5'UTR	   7
# 8	YKL109W	  5'UTR	  239

###################################compare A site mapping#########################

#####5UTR:
hdf5file <- rhdf5::H5Fopen(hd_file_none)
meta_5genes_none <- UTR5_table(genes)
##just do the plotting here

hdf5file <- rhdf5::H5Fopen(hd_file_CHX)
meta_5genes_CHX <- UTR5_table(genes)

hdf5file <- rhdf5::H5Fopen(hd_file_3AT)
meta_5genes_3AT <- UTR5_table(genes)

#for WT_none
meta_5genes_plot_none <- plotting_meta_analysis(meta_5genes_none)
ggsave("5UTR_notAmapped_none", device = "jpg")

#for WT_CHX
meta_5genes_plot_CHX <- plotting_meta_analysis(meta_5genes_CHX)
ggsave("5UTR_notAmapped_CHX", device = "jpg")

#for WT_3AT
meta_5genes_plot_3AT <- plotting_meta_analysis(meta_5genes_3AT)
ggsave("5UTR_notAmapped_3AT", device = "jpg")



###################################G-Sc_2014#######################################

#####1) A-site mapped:
WT_none <- All_genesAmapped(gene = genes, dataset = dataset, hdf5file = hdf5file_none, gffdf = gff_df, min_read_length = 10)

none_plot <- plotting_5UTR_with_title(WT_none)
ggsave("A-site mapped none", device = "jpg")

WT_CHX <- All_genesAmapped(gene = genes, dataset = dataset, hdf5file = hdf5file_CHX, gffdf = gff_df, min_read_length = 10)

CHX_plot <- plotting_5UTR_with_title(WT_CHX)
ggsave("A-site mapped CHX", device = "jpg")

WT_3AT <- All_genesAmapped(gene = genes, dataset = dataset, hdf5file = hdf5file_3AT, gffdf = gff_df, min_read_length = 10)

AT3_plot <- plotting_5UTR_with_title(WT_3AT)
ggsave("A-site mapped 3AT", device = "jpg")


###############plotting 3 conditions at the same time
none <- All_genesAmapped(gene = gene_names, dataset = dataset, hdf5file = hdf5file_none, gffdf = gff_df, min_read_length = 10)

CHX <- All_genesAmapped(gene = gene_names, dataset = dataset, hdf5file = hdf5file_CHX, gffdf = gff_df, min_read_length = 10)

AT3 <- All_genesAmapped(gene = gene_names, dataset = dataset, hdf5file = hdf5file_3AT, gffdf = gff_df, min_read_length = 10)

plot_none <- plotting_5UTR(none)
plot_CHX <- plotting_5UTR(CHX)
plot_AT3 <- plotting_5UTR(AT3)

ggarrange(plot_none,plot_CHX,plot_AT3, labels = c("None", "CHX", "AT3"), ncol = 1, nrow = 3)
ggsave("G-Sc_2014 3 conditions", device = "jpg")

# single gene

YOR303W_CHX <-A_mapped_genes("YOR303W", dataset = dataset_G2014, hdf5file = hdf5file_CHX, gffdf = gff_df, min_read_length = 10)

plotting_multiple(YOR303W_CHX, "YOR303W")
ggsave("YOR303W_CHX", device = "jpg")

something_none <-GetGeneCodonPosReads1dsnap(gene = "YEL009C", dataset = dataset, hdf5file = hdf5file_none, gffdf = gff_df, nnt_gene = nnt_gene, min_read_length = 10, asite_disp_length = asite_disp_length)

something_CHX <-GetGeneCodonPosReads1dsnap(gene = "YEL009C", dataset = dataset, hdf5file = hdf5file_CHX, gffdf = gff_df, nnt_gene = nnt_gene, min_read_length = 10, asite_disp_length = asite_disp_length)

something_3AT <-GetGeneCodonPosReads1dsnap(gene = "YEL009C", dataset = dataset, hdf5file = hdf5file_3AT, gffdf = gff_df, nnt_gene = nnt_gene, min_read_length = 10, asite_disp_length = asite_disp_length)

something_n_plot <- plotting_5UTR(something_none)
something_C_plot <- plotting_5UTR(something_CHX)
something_3_plot <- plotting_5UTR(something_3AT)

ggarrange(something_n_plot,something_C_plot,something_3_plot, labels = c("None", "CHX", "AT3"), ncol = 1, nrow = 3)
ggsave("YEL009C 3 conditions", device = "jpg")

## none
NONE_sev <- A_mapped_genes(genes, dataset = dataset, hdf5file = hdf5file_none, gffdf = gff_df, min_read_length = 10)


NONE_sev_5UTR <- CanVsNon(NONE_sev)

NONE_sev_5UTR_scatter <- ggplot(NONE_sev_5UTR, aes(x = Counts_AUG, y = Counts_UTR5)) +
  geom_point() + 
  geom_text(aes(label = Gene))
ggsave("NONE_sev_5UTR_scatter", device = "jpg")

NONE_regions <- all_regions_together(NONE_sev, gene_names = genes)


## CHX
CHX_sev <- A_mapped_genes(genes, dataset = dataset, hdf5file = hdf5file_CHX, gffdf = gff_df, min_read_length = 10)

CHX_sev_5UTR <- CanVsNon(CHX_sev)

CHX_sev_5UTR_scatter <- ggplot(CHX_sev_5UTR, aes(x = Counts_AUG, y = Counts_UTR5)) +
  geom_point() + 
  geom_text(aes(label = Gene))
ggsave("CHX_sev_5UTR_scatter", device = "jpg")


CHX_regions <- all_regions_together(CHX_sev, gene_names = genes)

## 3AT
AT3_sev <- A_mapped_genes(genes, dataset = dataset, hdf5file = hdf5file_3AT, gffdf = gff_df, min_read_length = 10)

AT3_sev_5UTR <- CanVsNon(AT3_sev)

AT3_sev_5UTR_scatter <- ggplot(AT3_sev_5UTR, aes(x = Counts_AUG, y = Counts_UTR5)) +
  geom_point() + 
  geom_text(aes(label = Gene))
ggsave("AT3_sev_5UTR_scatter", device = "jpg")
               
AT3_regions <- all_regions_together(AT3_sev, gene_names = genes)

###### comparisons

compared_Guydosh <- UTR5_3_conditions_all(NONE_sev_5UTR, CHX_sev_5UTR, AT3_sev_5UTR)
ggsave("compared_Guydosh", device = "jpg")
#sprawdzone ze kolejosc jest w porzadku 

### 5UTR comparison for one gene genes (with META_I_1_sev_5UTR result)
compared_UTR_Guydosh_YML065WC <- compared_UTR_single(NONE_sev_5UTR, CHX_sev_5UTR, AT3_sev_5UTR, "YML065W") 
ggsave("compared_UTR_YML065WC", device = "jpg")


