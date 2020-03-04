###Set up:

nnt_gene<- 50
setwd("/Users/Ania/Desktop/Szkoła/4th year/Dissertation/gene_graphs/")

dataset <- "G-Sc_2014"
hd_file <- "G-Sc_2014/output/WTnone/WTnone.h5"

# prepare files, opens hdf5 file connection
hdf5file <- rhdf5::H5Fopen(hd_file) # filehandle for the h5 file

gene_names <- rhdf5::h5ls(hdf5file, recursive = 1)$name
# next step is to test with these 5 genes only instead of all genes
test_orfs<- c("YCR012W","YEL009C","YOR303W","YOL130W","YGR094W")

orf_gff_file <- "G-Sc_2014/input/yeast_CDS_w_250utrs.gff3"

gff_df <- readGFFAsDf(orf_gff_file)

####Functions

readGFFAsDf <- purrr::compose(
  rtracklayer::readGFFAsGRanges,
  data.frame, 
  as_tibble,
  .dir = "forward" # functions called from left to right
)

GetGeneDatamatrix <- function(gene, dataset, hdf5file) {
  hdf5file %>%
    rhdf5::H5Dopen(
      name = paste0("/", gene, "/", dataset, "/reads/data")
    ) %>%
    rhdf5::H5Dread() %>%
    return()
}

#two functions because GetGeneDatamatrix doesnt accept -ve values 
Get_Get5UTRstart <- function(name, gffdf, ftype="UTR5", fstrand="+") {
  gffdf %>% 
    dplyr::filter(type==ftype, Name == name, strand == fstrand) %>% 
    dplyr::pull(start)
}

Tidy_Get5UTRstart <- function(name, gffdf, ftype="UTR5", fstrand="+") {
  gffdf %>% 
    dplyr::filter(type==ftype, Name == name, strand == fstrand) %>% 
    dplyr::pull(UTR5start)
}

UTR5_length_table <- function(gffdf, ftype="UTR5", fstrand="+") {
  gffdf %>% 
    filter(type==ftype, strand == fstrand)} 

#### Modification of the table for the GetGeneDataMatrix(gff_df) and TidyMatrix function (gff_df2)

table_modifications<- function(gffdf, ftype = "UTR5", fstrand = "+", End = "end", Start = "start") {
gffdf %>%
    filter(type==ftype, strand == fstrand) %>%
  mutate(UTR5start = gffdf$End *-1, UTR5L = seq(gffdf$Start:gffdf$End), UTR5_end = gffdf$Start-2 ) %>%
  select()
} ##okay it's clearly not liking this function bc of the End and Start. I want to avoid gff_df$. but once it works out it can replace the commands below

gff_df <- UTR5_length_table(gff_df)

#this is separate for tidydatamatrix purpouses
gff_df2 <- mutate(gff_df, UTR5L=gff_df$end-gff_df$start)
gff_df2 <- mutate(gff_df, UTR5start = -(gff_df$end), UTR5_end = gff_df$start-2)
gff_df2 <- gff_df2[ -c(2:4) ] #for the moment let's say it's correct

##back to the functions 
UTR5_length <- function(name, gffdf, ftype="UTR5", fstrand="+") {
  gffdf %>% 
    dplyr::filter(type==ftype, Name == name, strand == fstrand) %>% 
    dplyr::pull(UTR5L)
}

GetGeneDatamatrix5start <- function(gene, dataset, hdf5file, 
                                    Get_Get5UTRstart, UTR5full, nnt_gene) {
  data_mat_all <- GetGeneDatamatrix(gene, dataset, hdf5file)
  #I'm changing n_left5 to just 1st nt of 5'UTR
  n_left5 <- Get_Get5UTRstart # column to start from (5'end)
  zeropad5_mat <- matrix(0, nrow = nrow(data_mat_all), ncol = 0)
  #or ncol = (n_buffer - posn_5start + 1)
  n_right3 <- UTR5full + nnt_gene - 1 # column to end with (3'end) 
  data_mat_5start <- data_mat_all[, n_left5:n_right3]
  return(cbind(zeropad5_mat, data_mat_5start))
}


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

#all_genes_all_info <- as_tibble(interesting ) #I don't really get why that doesnt work, the only difference is that I am not specifying which gene even though it only contains 1 matrix. Am I opening up the matrix somehow with the second function? okay i think i know,the gene must be selected or otherwise it will turn all the lists into one tibble which isnt what we want 

#here i need to look into a function that separates separate matrix for each gene? or that at least treates each matrix separately 

wtf<-as_tibble(interesting$YCR012W, .name_repair ="minimal")

UTR5_YCR012W<-Tidy_Get5UTRstart(test_orfs[1], gff_df2)

positions <- UTR5_YCR012W : (UTR5_YCR012W +ncol(wtf)-1)

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


##Selecting 5'UTR start position required for plotting
#the function would identify the first position which identifies read counts and from that (or    this position+ for instance 10 nt just for visualization)

