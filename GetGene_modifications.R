rm(list=ls())
if (!require(xlsx)){
  install.packages("xlsx")
  library(xlsx)
}

#Planning:

#basically the GetGeneDataMatrix seems fine:
GetGeneDatamatrix <- function(gene, dataset, hdf5file) {
  hdf5file %>%
    rhdf5::H5Dopen(
      name = paste0("/", gene, "/", dataset, "/reads/data")
    ) %>%
    rhdf5::H5Dread() %>%
    return()
}

#it just opens up the file and gives us datasets within.
#but then we have to limit these to the ones we want
  #we want: 5'UTR and for the start the entirety of the CDS, it's just easier that way

# > gff_df
# # A tibble: 17,436 x 10
# seqnames  start   end width strand source      type  score phase Name     
# <fct>     <int> <int> <int> <fct>  <fct>       <fct> <dbl> <int> <chr>    
# 1 YAL068C       1   250   250 +      rtracklayer UTR5     NA    NA YAL068C  
# 2 YAL068C     251   613   363 +      rtracklayer CDS      NA    NA YAL068C  
# 3 YAL068C     614   863   250 +      rtracklayer UTR3     NA    NA YAL068C

Get5UTRstart <- function(name, gffdf, ftype="UTR5", fstrand="+") {
  gffdf %>% 
    dplyr::filter(type==ftype, Name == name, strand == fstrand) %>% 
    dplyr::pull(start)
}

Get5UTRstart(name = test_orfs[1],gffdf = UTR5_length_table_output)

# Changed ftype="CDS" to ftype="UTR5"

#checked that Get5UTRstart works

##5UTR length--> required to calculate n_right3 (I think at least)

UTR5_length_table <- function(gffdf) {
  gffdf %>%
  mutate(UTR5L=gff_df$end-gff_df$start)
}

#instead of creating a new table, just reassign it as gff_df==that's gonna make it much simpler
gff_df<-UTR5_length_table(gff_df)

UTR5_length <- function(name, gffdf, ftype="UTR5", fstrand="+") {
  gffdf %>% 
    dplyr::filter(type==ftype, Name == name, strand == fstrand) %>% 
    dplyr::pull(UTR5L)
}

#UTR5_length("YCR012W",UTR5_length_table)
  #output: 249 "CORRECT"

#now the issue us to use it in the next command-we dont really want the 5'UTR to be that long if it's not necessary 
  #it would have to be individual to different sequences: 
    #maybe the function could somehow count the threshold or something? 
    #or the first column could start from the first footprint found<--this one sounds good 

#GetGeneDataMatrix5Start: 
#so this is supposed to give us the table from 5'UTR to the end of the CDS

GetGeneDatamatrix5start <- function(gene, dataset, hdf5file, 
                                    Get5UTRstart, UTR5full, nnt_gene) {
  data_mat_all <- GetGeneDatamatrix(gene, dataset, hdf5file)
  #I'm changing the n_left to just 1st nt of 5'UTR
      n_left5 <- Get5UTRstart # column to start from (5'end)
    zeropad5_mat <- matrix(0, nrow = nrow(data_mat_all), ncol = 0)
   
  n_right3 <- UTR5full + nnt_gene - 1 # column to end with (3'end) 
  data_mat_5start <- data_mat_all[, n_left5:n_right3]
  return(cbind(zeropad5_mat, data_mat_5start))
}

interesting <-
  lapply(test_orfs[1],
         function(gene) 
           GetGeneDatamatrix5start(gene,
                                   dataset,
                                   hdf5file,
                                   Get5UTRstart(gene, gffdf = gff_df),
                                   UTR5full = UTR5_length(gene, gffdf = gff_df), 
                                   nnt_gene = 50)
  ) %>%
  Reduce("+", .) %>% # sums the list of data matrices
  TidyDatamatrix(startpos = -250 + 1, startlen = 10)


#maybe I can assume that the reduce bit is fine
#and then the TidyDataMatrix

#why is the startpos and startlen = 1?
TidyDatamatrix <- function(data_mat, startpos = 1, startlen = 1) {
  # CHECK startpos/off-by-one
  positions <- startpos:(startpos + ncol(data_mat) - 1)
  readlengths <- startlen:(startlen + nrow(data_mat) - 1)
  data_mat %>%
    set_colnames(positions) %>%
    as_tibble() %>%
    mutate(ReadLen = readlengths) %>%
    gather(-ReadLen, key = "Pos", value = "Counts", convert = FALSE) %>%
    mutate(Pos = as.integer(Pos), Counts = as.integer(Counts))
}

### the functio

##Selecting 5'UTR start position required for plotting
  #the function would identify the first position which identifies read counts and from that (or    this position+ for instance 10 nt just for visualization)


####Plotting

interesting1<-interesting %>%
  group_by(Pos) %>%
  summarise(Counts=sum(Counts)
  )

interesting1<- arrange(interesting1,interesting1$Pos)

#Plotting 

interesting1_plot<-ggplot(
  data=interesting1,) +
  geom_histogram(aes(x=Pos, y=Counts), stat="identity") + 
  theme_light() +
  labs(y= "Read count", x = "Position (mRNA)")

ggsave("Ania",plot = interesting1_plot, device = "png")



### Potential functions?

Find() #returns the first element which matches the predicate (or the last element if right = TRUE).
Position() #returns the position of the first element that matches the predicate (or the last element if right = TRUE).


integrate() #finds the area under the curve defined by f()
uniroot() #finds where f() hits zero
optimise() #finds the location of lowest (or highest) value of f()


### Playing
gene <- as.list(test_orfs[1])
gene %>%
GetGeneDatamatrix5start(gene,
                        dataset,
                        hdf5file,
                        Get5UTRstart= Get5UTRstart(gene, gffdf = gff_df),
                        UTR5full = UTR5_length(gene, gffdf = gff_df), 
                        nnt_gene = 50)
traceback()




 %>%
  Reduce("+", .) %>% # sums the list of data matrices
  TidyDatamatrix(startpos = -250 + 1, startlen = 10)

startpos:(startpos + ncol(data_mat) - 1)

play

#we don't want the reduce function at all-bc it compresses different genes intothe same positions
#instead I need to find a way to modify the TidyDataMatrix without the reduce function

  

TidyDatamatrix <- function(data_mat, startpos = 1, startlen = 1) {
  # CHECK startpos/off-by-one
  positions <- startpos:(startpos + ncol(data_mat) - 1)
  readlengths <- startlen:(startlen + nrow(data_mat) - 1)
  data_mat %>%
    set_colnames(positions) %>%
    as_tibble() %>%
    mutate(ReadLen = readlengths) %>%
    gather(-ReadLen, key = "Pos", value = "Counts", convert = FALSE) %>%
    mutate(Pos = as.integer(Pos), Counts = as.integer(Counts))
}