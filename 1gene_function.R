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
#WT_none <- "G-Sc_2014/output/WTnone/WTnone.h5"
WT_3AT <- "G-Sc_2014/output/WT3AT/WT3AT.h5"
#WT_CHX <- "G-Sc_2014/output/WTCHX/WTCHX.h5"

# prepare files, opens hdf5 file connection
dataset <- "G-Sc_2014"
hd_file <- WT_3AT
hdf5file <- rhdf5::H5Fopen(hd_file) # filehandle for the h5 file

#Initial set ups
test_orfs <- c("YCR012W","YEL009C","YOR303W","YOL130W","YGR094W")
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
    # uAUG_efficiency %>%  optional to put it here
  return()
}


#the final function i'd use :) iteration of UTR5_table 
final_function_table <- function(x){
 purrr::map(x, UTR5_table)  ### as_tibble(.name_repair = "unique)
}

output_orfs<-final_function_table(test_orfs) #it works just fine 


############################################################################################                                               ##plotting##

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
#
#                                     ##multiple plots##
#
# ###ok so that's the final function for plotting all graphs from a list of tibbles (output from )
#     #still dont know how to name each graph and give a title to the plot overall
plotting_multiple <- function(input_data) {
  purrr :: map(input_data, plotting_5UTR) %>%
    ggarrange(plotlist = .) %>%  ##, common.legend = (maybe to establish a name)
    return()
    #ggsave(file.path(paste0("Plot of ", gene)), device = "jpg")
}
#
# xxx<-plotting_multiple(output_orfs)
# #issue with the name now..same for all, i think it means we need to name them earlier on so we know what's what. Why does it only use the YCR012W?
#
# ggsave("multiple", device = "jpg")

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

###################################################################################
                              #A site displacement slot#

###
CalcAsiteFixedOneLength <- function(reads_pos_length, min_read_length,
                                    read_length, asite_disp) {
  # Calculate read A-site using a fixed displacement for a single read length
  length_row_choose <- read_length - min_read_length + 1
  reads_pos_length[length_row_choose, ] %>%
    dplyr::lag(n = asite_disp, default = 0)
}

something <- CalcAsiteFixedOneLength(reads_pos_length = GetGeneDatamatrix(gene = test_orfs[1],dataset = dataset,hdf5file = hdf5file), min_read_length = 10, read_length = 28, asite_disp = 15)
# > str(something)
# num [1:1751] 0 0 0 0 0 0 0 0 0 0 ...

# > something
# [1]    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
# [20]    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
# [39]    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
# [58]    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
# [77]    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
# [96]    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
# [115]    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
# [134]    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
# [153]    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
# [172]    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    1    0    0
# [191]    0    0    0    0    0    0    1    0    0    1    0    0    2    0    0    0    3    0    0
# [210]    1    0    9    1    0    1    0    0    4    5    0    0    0    0    0    0    0    0    0
# [229]    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
# [248]    0    0    0    0    0    0   80    8    0    6    2    0    6    1    0    9    2    0   59
# [267]   18    0    6    6    3  129    6    1  278   20    0   40   17    3   42    3    0    1    1
# [286]    1   19    0    0    0    4    9   36   14    0   42    6   33   25    2    0   12    1    0
# [305]   22    1    0    1    0    0  162    0    3   71   51    0   73   53   18   55   25    1   37
# [324]    3    0  502   65   10   83   54    6   10   20    7   19    2    0   53    1    0   11    0
# [343]    1   27   53    1  279   30    1  283   13    0    0   23    0    0    3    0   10    0    1
# [362]    7    3    0   62   17    2  120   62    4  107   21   21   10    5    2   40  127    0  101
# [381]   12    5   91   53    1   33   34    2    8   11   12  139   28   16   13    9    0   87   32
# [400]    0  114   18    2  352   30   13 2532   28   13   86  137    2  273   70   38   61   23    1

###
CalcAsiteFixed <- function(reads_pos_length, min_read_length,
                           asite_disp_length = data.frame(
                             read_length = c(28, 29, 30),
                             asite_disp = c(15, 15, 15)
                           ),
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

CalcAsiteFixed(reads_pos_length =GetGeneDatamatrix(gene = test_orfs[1],dataset = dataset,hdf5file = hdf5file), min_read_length = 10 , asite_disp_length =data.frame(
  read_length = c(28, 29, 30),
  asite_disp = c(15, 15, 15)
) , colsum_out = TRUE)

# [1]    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
# [17]    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
# [33]    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
# [49]    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
# [65]    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
# [81]    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
# [97]    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
# [113]    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
# [129]    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
# [145]    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
# [161]    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
# [177]    0    0    0    0    0    0    0    0    0    0    0    1    0    1    0    0
# [193]    0    0    0    0    1    0    0    1    0    8    2    0    0    1    3    0
# [209]    7    2    5    9    1    1    2    0    2   78    9    0    0    3    0    0
# [225]    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0
# [241]    0    0    0    1    0    0    0    0    0    0    0    0    0  134    8    1
# [257]   17    2    0   18    1    1   10    2    3   59   20    0    8    7   89  130
# [273]    6    2  386   21   20   42   17    4   44    3    0    6    1    1   27    0
# [289]    5    3    4   18  142   14    5   72    8  188   29    8    0   28    1    1
# [305]   28    1    0    7    2   31  216    1    3  105   54    8  182   57  204   79
# [321]   29    6  245    5    1  535   68   93 1192   76   10   80   29   23   21    3
# [337]    1   57    2    1   15    0   23  407   78    8  312   32    3  485   15    0
# [353]   22   23    0    0    4    2   12    0    1   14    3    4   69   19    2  123
# [369]   64   13  109   25   79   21    5    6   82  157    0  179   14   21  249   61
# [385]   21   36   36   11    8   21  412  144   37   99   19   11    0  103   32    1
# [401]  136   20   11  412   65   37 3106   33   48  904  144    8  350  128  293  102
# [417]   30   15  224   36    1   10    0    0    3    1    1   31   99    4  198    2
# [433]   13  205   10  132 2466  710   28   93   98    4   26   17    2  100   41   45
# [449]   19    5   11   62    7    1  448   22   12   81   29  197   53   48    3   61
# [465]    5   60  304   41    4   25   15   54    5   44    2  137    6   29  124    7
# [481]   10   66   29  331  122    3    4   21    2    0  104    1  364   38   17    3
# [497]   58    4    2    0    0    2   77   11   61   82   19   77   71   19   11   78
# [513]   32    3  117    4   30  530  148    3  115   21    1   31    2    1   29    4
# [529]   14   64   18   89   72   28    5  115   11  205  568   28    9   50   21    3
# [545]   67    7    0   95    3    0   81    2    2  195    9    4   82    4   90   21
# [561]   33    4  189   34   60   81   15    1  174   23    0   25    1    1    8    2
# [577]    1   32    5    1   34    7   25   73    5    1   78   11   15   50   20   15
# [593]  132   17  184   15   42    3   16   36    1   16   24   43   11    3    5  142
# [609]   38   22  135   30   11   18    6    1  105    4   31  549  245   21  107   15
# [625]    6   67    7    0   17   10   83   15    0   20  213    5   20  486   28   23
# [641]  290   60    7   32    7    0   30    4    5   37   17   18  121   14    3  235
# [657]   22   12    3    5    2   49    4    2  180   35  224   47   50   53  193  124
# [673]   15  116    9    2   33    2    1   40    7    2   10    7   61   67    2    5
# [689]   37    4    4  122   26    6   89   11  846  115  170   10   63   14   14   41
# [705]   19   22  549   25   17  105    7    5   57    2    2   16    4    1   44    8
# [721]    5  206   22   10   35    9    9    2    1    6  444   25  218   30   37   68
# [737]  137   45   38   39    7    3  372   15   19  545   12    5  206   79 3445 2647
# [753]  146   40  445   91    3   14    9  202  689  119   19   72  137  164  231   22
# [769]    4  307   74   20  132   46   86   42   17    7   67    4   15   50    1   57
# [785]   14    5    0  112   36    7  116   19   59    8    6    1   17    1    1   47
# [801]    3    5   63   17    7  171  118    6  231   93   62  164   44   29   24   22
# [817]    0   16    7  195  111    3    3   40    9    9  424   16   11  216   82   61
# [833]  462    7    5   69   16    1  105    0    0  272    2    3  141   19    9  189
# [849]    6    8  185   10    4    9    0    2  141    7   31   49    9    2   70    8
# [865]   28  109    5    1   39   18    3   37   15    2    9    0    0   11    3    0
# [881]   40    0    1  161    4    1   46    2    0   12    0    0   12    0    0    7
# [897]    6   31   25  190   50   86   30   12  262    2   16  523   82   60   19   37
# [913]   31   47    4    0  149    9   24   95   71  183  104  179    2  226    7    7
# [929]  284  193   82   17   28    2  135    0    2  173    6    4  177    8    2   43
# [945]    3    3    2    0    2  185    4    1   27   15    1   44    1   14   92    9
# [961]   11  228   18    2   40   24   66   21   34  187  184  210   11   26   48   21
# [977]  180   33    6   32   14   14   15    1    0   25    3    2  159    1    1   13
# [993]    0    0   53    2    1  110   21    4
# [ reached getOption("max.print") -- omitted 751 entries ]

###
SnapToCodon <- function(x, left, right, snapdisp=0L) {
  # snap nucleotide-aligned reads to codon position
  #   x:     vector
  #   left:  integer for starting position, frame 0
  #   right: integer for ending position
  #   snapdisp: integer any additional displacement in the snapping
  RcppRoll::roll_suml(x[(left:right) + snapdisp], n=3L, by=3L, fill = NULL)
}

left = getCDS5start(test_orfs[1], gff_df) #finds left and right value for each gene 
right = GetCDS3end(test_orfs[1], gff_df)

getCDS5start <- function(name, gffdf, ftype="CDS", fstrand="+") {
  gffdf %>% 
    filter(type==ftype, Name == name, strand == fstrand) %>% 
    pull(start) %>% 
    min 
}

GetCDS3end <- function(name, gffdf, ftype="CDS", fstrand="+") {
  gffdf %>% 
    filter(type==ftype, Name == name, strand == fstrand) %>% 
    pull(end) %>% 
    max 
}


###
GetGeneCodonPosReads1dsnap <- function(gene, dataset, hdf5file, left, right, 
                                       min_read_length, 
                                       asite_disp_length = data.frame(
                                         read_length = c(28, 29, 30),
                                         asite_disp = c(15, 15, 15)
                                       ), 
                                       snapdisp=0L) {
  reads_pos_length <- GetGeneDatamatrix(gene, dataset, hdf5file) # Get the matrix of read counts
  reads_asitepos <- CalcAsiteFixed(
    reads_pos_length, min_read_length,
    asite_disp_length
  )
  SnapToCodon(reads_asitepos,left,right,snapdisp)
}

final_A_mapped <- GetGeneCodonPosReads1dsnap(gene = test_orfs[1], dataset = dataset ,hdf5file = hdf5file ,left =left , right = right , min_read_length = 10 , asite_disp_length = data.frame(
  read_length = c(28, 29, 30),
  asite_disp = c(15, 15, 15)), snapdisp = 0L)

# #[1]    0  143   19   20   15   79  104  138  427   63   47    8   32   25  161  268   37   30   29
# [20]   40  220  167  443  114  251  696 1278  132   25   60   38  493  347  500   45    6   13   21
# [39]   90  200  213   32  239  214  331   83  441  280   30  136  167  514 3187 1056  771  147  261
# [58]   10    5  134  213  347 3204  195   45  186   35   70  482  307  104  126  349   94   51  172
# [77]  141  426  129   23  469   58   64    2  149  178  101  113  151  681  137   34   47  171  105
# [96]  331  605   74   74   98   85  208  176   58  283   97  197   27   11   38   66   79  104   85
# [115]  333   60   53   83   19  202  176   25  140  815  128   74  110   35  238  537  357   39   39
# [134]   72  138  269   10   55  439  150  332  127   36   49   78   74   45  154  946  295   91   82
# [153]  591  117   61   21   57  238   53    9  687  135  220   49  406  562 3730 2833  539  225  827
# [172]  373  257  401  264   66   86  108   19  155  194   15   19   55   87  295  386  237   46  218
# [191]  117   58  451  359  474   86  105  277  169  203  199   11  179   60  106  115   60   54    9
# [210]   14   41  166   48   12   12   44  265  128  280  665   87   51  182  349  285  240  559   47
# [229]  137  183  187   49    4  190   43   59  112  248  130  242  405   95  219   60   16   30  161
# [248]   13   56  135    4  150  113  353  312  378  145   34   31   79   35   53  217  284  244  159
# [267]   67  236  237  350  539  105  361  635  573  254  519  160  187  326  222   26   17  100  582
# [286]   92   91   98  212  478   41  175   96  538   51  626  245  255  148  163   27  173   49  154
# [305]  320  135   50  219  584  200  313  332  109   72   13   18  148  440  315   33   24   75  205
# [324]  185  459  219  177   48  551  307  163  123  278  111   62  110  115  192 1468  463  602  443
# [343]  147   64  366   99   69  213  121  154  203  318  193   33   39   51   92   36  134   45  162
# [362]  312  243   58   39   45   40   30   13   46   21  189   16  173  385  302  110  248  523   83
# [381]  250  297  418  561  375  429  318  202 5581  644   31   12   29   36  109   62   69   92  191
# [400]   89  148  317   36   69  107  136  180  109   34   78  172  478  293   41   25   12   24
> 



##########################################################################################
                                      ##Zooming in ## 

# test object with 1 gene
#tibble1 <- output_orfs[[1]]
# > str(tibble1)
# Classes ‘tbl_df’, ‘tbl’ and 'data.frame':	300 obs. of  2 variables:
#   $ Pos   : int  -250 -249 -248 -247 -246 -245 -244 -243 -242 -241 ...
# $ Counts: int  0 0 0 0 0 0 0 0 0 0 ...

# tibble5 <- output_orfs[[5]] i've created more tables to play with 


##### Function for the upstream codons 

# finding_uAUG_beginning <- function(datatibble) {
#   datatibble %>% 
#     dplyr :: filter(Counts > 0) %>%
#     min() %>%
#     return()
# }

###i was playing with tibble5 (instead of tbl1) bc it has much lower count frequency, i encountered the following issue: 
#for finding_uAUG_beginning(tibble5) if filter(Counts > 4) bc the max count value of tbl5 is 4 so that's a good question how i can tackle that. 
# Error in FUN(X[[i]], ...) : 
#   only defined on a data frame with all numeric variables 

# to run: 
#finding_uAUG_beginning(tibble1)
# example run, set to uAUG_beginning:
#uAUG_beginning <- finding_uAUG_beginning(tibble1)

#finding_uAUG_ending <- function(uAUG_start) {
#   uAUG_start + temporary_length %>%
#     return()
# }

# to run: 
#finding_uAUG_ending(uAUG_beginning)
# example run, set to uAUG_ending:
#uAUG_ending <- finding_uAUG_ending(uAUG_beginning)

sum_uAUG <- function(datatibble, uAUG_start, uAUG_end){
  datatibble %>%
    filter(Pos >= uAUG_start, Pos <= uAUG_end) %>%
    summarise(sum_read_counts_uAUG=sum(Counts)) %>%
    return()
}

# to run: 
#sum_uAUG(tibble1, uAUG_beginning, uAUG_ending)
# example run set to sum_uAUG_result
#sum_uAUG_result <- sum_uAUG(tibble1, uAUG_beginning, uAUG_ending)


####AUG
finding_AUG_beginning <- function(datatibble) {
  datatibble %>% 
    dplyr::select(Pos) %>%
    filter(Pos == -10) %>%  #I changed the value to 6 so altogether now it's 12nt region
    pull 
}

# to run: 
#finding_AUG_beginning(tibble1)
# example run, set to uAUG_beginning:
#AUG_beginning <- finding_AUG_beginning(tibble1)

finding_AUG_ending <- function(AUG_start) {-
    AUG_start + temporary_length %>%
    return()
}

# to run: 
#finding_AUG_ending(AUG_beginning)
# example run, set to uAUG_ending:
#AUG_ending <- finding_AUG_ending(AUG_beginning)

sum_AUG <- function(datatibble, AUG_start, AUG_end){
  datatibble %>%
    filter(Pos >= AUG_start, Pos <= AUG_end) %>%
    summarise(sum_read_counts_AUG=sum(Counts)) %>%
    return()
}
# to run: 
#sum_AUG(tibble1, AUG_beginning, AUG_ending )
# example run, set to uAUG_ending:
#sum_AUG_result <-sum_AUG(tibble1, AUG_beginning, AUG_ending )

##efficiency
upstream_efficiency <- function(uAUG_sum_result, AUG_sum_result) {
  (uAUG_sum_result/AUG_sum_result *100) %>% 
    return()
}
  
#upstream_efficiency(sum_uAUG_result/sum_AUG_result)

##final function 

uAUG_efficiency <- function(datatibble, uAUG_start,uAUG_end) {
#finding uAUG and AUG endpoints
  AUG_start <- -10
  AUG_end <- AUG_start + temporary_length
  # AUG_start <-finding_AUG_beginning(datatibble)
  # AUG_end <- finding_AUG_ending(AUG_start)
  # uAUG_start <-finding_uAUG_beginning(datatibble)
  # uAUG_end <- finding_uAUG_ending(uAUG_start)

  #summing uAUG and AUG region
  final_sum_uAUG <- sum_uAUG(datatibble, uAUG_start, uAUG_end )
  final_sum_AUG <- sum_AUG(datatibble, AUG_start, AUG_end )
  
  #calculating the efficiency
  final_result_of_efficiency <- upstream_efficiency(final_sum_uAUG, final_sum_AUG) 
  
  paste0("Efficiency of upstream translation initiation is ", final_result_of_efficiency, " % relative to the AUG start site with " , final_sum_uAUG , " read counts in the UTR and ", final_sum_AUG , " read counts around the AUG" )
}

##
uAUG_efficiency_many_genes <- function(datatibble, uAUG_start, uAUG_end) {
  map(datatibble, uAUG_efficiency, uAUG_start, uAUG_end)
}
## example:
uAUG_efficiency_many_genes(output_orfs,uAUG_start = -Inf,uAUG_end = -11 ) #=> output as expected 


####################################################################################################
                                   ##Working space for tibble5##

##idea: find different positions at which read count > 0, calculate the efficiency of that region relative to the AUG site, (which must be the same for AUG and each non-AUG), start counting from counts > 0: if larger than some % threshold then calculate the efficiency, if smaller just skip, do that for each counts > 0 result, the output will be only done for sites with larger thresholds

##### Function for the AUG start codons needs to be changed
# tibble5 <- output_orfs[[5]]
# 
# finding_AUG_beginning <- function(datatibble) {
#   datatibble %>%
#     dplyr::select(Pos) %>%
#     filter(Pos == -10) %>%  #I changed the value to 6 so altogether now it's 12nt region
#     pull
# }

#so i'm first calculating the efficiency of AUG initiation
#start of the search, i could be just looking at a specific region flanking the AUG site, maybe it could be maybe 24 NTs?

#step: find each counts > 0 and calculate the % for it

# you're getting a table with only counts > 0 values,
# finding_uAUG_beginning2 <- function(datatibble) {
#   datatibble %>%
#     dplyr :: filter(Counts > 0) %>%
# 
# i could create a function that takes the uAUG positions and creates a vector of every 20th position, then we select only these raws
#   > seq(from=0,to=100, by=20)
# [1]   0  20  40  60  80 100
# again it's not exactly the same, we're taking basically omitting every 20th value from the table
# another option
# m[seq(1, length(m), 20)] %>%good idea but not that exactly
#     return()
# 
# }
  
# finding_uAUG_beginning2_single <- function(datatibble) {
#   datatibble %>% 
#     dplyr :: filter(Counts > 0) %>%
#     min
# }
# uAUG_first_value <-finding_uAUG_beginning2_single(tibble5)
# # [1] -31
# 
# seq(from =uAUG_first_value, to)

#find beginning for each count > 0 site
# AUG_beginning_tibble5 <- finding_uAUG_beginning2(tibble5)
# > AUG_beginning_tibble5
#[1] -9 34 38


#find end for each count > 0 site
# AUG_end_tibble5<-finding_uAUG_ending(AUG_beginning_tibble5)
#[1]  21 -22 -26

#then calculate the percentage relative to AUG 


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
