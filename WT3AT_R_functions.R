dataset_WT3AT <- "G-Sc_2014"
hd_file_WT3AT <- "G-Sc_2014/output/WT3AT/WT3AT.h5"
hdf5file_WT3AT <- rhdf5::H5Fopen(hd_file) # filehandle for the h5 file

#Initial set ups
gene1 <- "YCR012W"
test_orfs <- c("YCR012W","YEL009C","YOR303W","YOL130W","YGR094W")
nnt_gene<- 50
startpos <-250
startlen <- 10
orf_gff_file <- "G-Sc_2014/input/yeast_CDS_w_250utrs.gff3"
gff_df <- readGFFAsDf(orf_gff_file)
WT3AT_data <-final_function(test_orfs,file = hdf)


### OMG IT WORKS!!! THIS IS AMAZING!!!