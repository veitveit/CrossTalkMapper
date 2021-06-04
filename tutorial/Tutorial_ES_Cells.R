##########################################################
##               CrossTalkMapper Tutorial               ##
## Analysis of mouse embryonic stem cells histones PTMs ##
##########################################################





###############=- TASK 1.1 -=##################
#-Packages installation and environment setup-#

##Vector of required packages' names
packages <-
  c(
    "broom",
    "tidyr",
    "gplots",
    "gtools",
    "ggplot2",
    "scales",
    "metR",
    "forcats",
    "ggrepel",
    "gridExtra",
    "dplyr",
    "paletteer"
  )
  
##If packages are missing run following lines
#update.packages(ask = FALSE)
#install.packages(setdiff(packages, rownames(installed.packages()))) 


##Load Packages
invisible(lapply(packages, library, character.only = TRUE))





#################=- TASK 1.2 -=###################
#-Set working directory and load CrossTalkMapper-#

##Define the working directory, by default is set to the location of this script:
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
##If it doesn't exist, makes a folder 'plots' in the working directory in which the generated plots will be stored,(note that a warning is given if the folder already exists)
dir.create(file.path(getwd(), "plots"))
##load PTM-CrosstalkMapper function (If "ptm-crosstalkmapper.R" is not located in the working directory please specify the path to the script)
source("../ctm-functions/ptm-crosstalkmapper.R") #TODO change path/script





##########=- TASK 2.1 -=###########
#-Import dataset and prepare data-#

## Import the dataset downloaded from CrosstalkDb and print the start of the dataframe to see if it requires formatting column organisation and name requirements ?
mytable <- read.csv("MouseStemCell.csv") #import the dataset downloaded from crosstalkDB, note that the file must be located in in current working directory define in Task 1.1
head(mytable)

## Here the structure and the column names of the dataframe "my_table" are adapted to fit the requirements of the function"prepPTMdata()"  
mytable$timepoint <- "None" # Create a column "timepoint" with "None" values 
mytable$tissue <- mytable$cell.type...tissue  # Duplicate the column "cell.type...tissue" in a column named "tissue"  
write.csv(mytable, "t.csv") # Save the formatted dataframe in csv as "t.csv"

## Prepare the formatted PTM data frame for the computation of PTM abundances and interplay score
data <- prepPTMdata(csv="t.csv",
                    histvars = TRUE, # histvars=TRUE is you wish to distinguish between H3.1 and H3.3. Else set to FALSE
                    avrepls = FALSE) # avrepls=TRUE is you wish to average replicates. Else set to FALSE





###############=- TASK 2.2 -=################
#-Compute PTM abundance and interplay score-#

ptm_ab <- calcPTMab(data) # This function calculates abundance of individual PTMs and co-occurence of PTM pairs as well as the interplay score for these PTM pair. 
ptm_ab$tissue <- factor( # Reorder the levels of the column "tissue" to display the WT condition first in plots 
    ptm_ab$tissue,
    levels = c("mESCs", "ESCs - Dnmt TKO", "ESCs - Ring1A1B-/-", "mESCs Suz12-/-")
  )




##############=- TASK 3.1 -=#################
#-Bar plot of single modification abundance-#

## A bar plot is created to visualize the distribution of individual PTMs and their abundance in the different conditions for the two histone variants H3.1 H3.3

## Extract individual PTM abundances from the PTM pairs abundances data frame and filter out the PTM with a low abundance  
plotdat <- unique(ptm_ab[, c("hist", "tissue", "timepoint", "repl", "mi", "pi")]) # Create a dataframe of single PTM abundances
filtered_ptms <- names(which(by(ptm_ab$pi, ptm_ab$mi, max) > 0.1)) # list of PTMs with abundance above threshold (0.1 in this case)
plotdat2 <- plotdat[plotdat$mi %in% filtered_ptms, ] # Filter dataframe based on the list 

## Create the bar plot from the filtered dataset 
p <- ggplot(plotdat2, aes(
  fill = tissue,
  x = reorder(mi,-pi), #Sort the bars in decreasing abundances 
  y = pi
)) +
  geom_bar(position = position_dodge(preserve = "single"), stat = "identity") +
  theme_minimal() +
  theme(text = element_text(size = 15),
        axis.text.x = element_text(
          angle = 90,
          vjust = 0.5,
          hjust = 1
        )) +
  coord_cartesian(clip = "off") +
  facet_wrap(~ hist) +
  xlab("Modification")

ggsave(filename = paste0("plots/ptm_distribution.pdf"), plot = p) #Save plot as PDF in the "plots" folder



##################################=- TASK 3.2 -=##########################################
#-heatmaps of Single PTM abundance, co-occurence, interplay score and peptide abundance-#


## This task presents how to create four different types of heatmaps 
## from a dataset of PTM abundances using the function heatmap_all() of CrossTalkMapper.

#-Heatmap of single PTMs abundance-#

## Prepare data
plotdat <- unique(ptm_ab[, c("hist", "tissue", "timepoint", "repl", "mi", "pi")]) # Extract abundance data for single PTMs 'mi'
flat_matrix <- # Convert the abundance dataframe to a flat_matrix
  reshape2::dcast(plotdat,
                  formula = mi ~ tissue + hist + timepoint + repl,
                  value.var = "pi")
## Plot data
heatmap_all(
  flat_matrix,
  showSidebar = "tissue", #  
  hscale = "none",
  title_of_plot = "H3 variants, single PTM frequency",
  label = "single_ptms",
  outdir = "plots/"
)


#-Heatmap of Co-occurring modification frequencies-#

## Prepare data
plotdat <- unique(ptm_ab[, c("hist", "tissue", "timepoint", "repl", "mi", "mj", "pi", "pij", "I")])
plotdat$mij <- paste(plotdat$mi, plotdat$mj, sep = "")
flat_matrix <-
  reshape2::dcast(plotdat,
                  formula = mij ~ tissue + hist + timepoint + repl,
                  value.var = "pij")
# Plot data
heatmap_all(
  flat_matrix,
  showSidebar = "tissue",
  hscale = "none",
  title_of_plot = "H3 variants, double PTM frequency",
  label = "double_ptms",
  outdir = "plots/"
)


#-Heatmap of Interplay scores-#

## Prepare data
ptm_ab_hist = ptm_ab[grepl("H3.1", ptm_ab$hist),]
plotdat <- unique(ptm_ab_hist[, c("hist", "tissue", "timepoint", "repl", "mi", "mj", "pi", "pij", "I")])
plotdat$mij <- paste(plotdat$mi, plotdat$mj, sep = "")
flat_matrix <-
  reshape2::dcast(plotdat,
                  formula = mij ~ tissue + hist + timepoint + repl,
                  value.var = "I")

## Plot data
heatmap_all(
  flat_matrix,
  showSidebar = "tissue",
  hscale = "none",
  title_of_plot = "H3 variants, interplay scores",
  label = "interplay_scores",
  outdir = "plots/"
)


#-Heatmap of peptides abundances#

## Prepare data (Note that in this case the unprocessed crosstalkDB dataset "mytable" is used, as the dataset "ptm_ab" does not contain quantification of peptides)
plotdat <-
  unique(mytable[, c("tissue", "modifications", "timepoint", "repl", "quantification")])
flat_matrix <-
  reshape2::dcast(
    plotdat,
    formula = modifications ~ tissue + timepoint + repl,
    value.var = "quantification",
    fun.aggregate = sum
  )

## Plot data
heatmap_all(
  flat_matrix,
  showSidebar = "tissue",
  hscale = "row",
  title_of_plot = "H3 peptides",
  label = "peptides",
  outdir = "plots/"
) 



###=- TASK 3.3 -=###
#-Crosstalk maps of K27 and K36 methylations-#

## Filter dataset to retain PTM pairs K27/K36 with any of the 3 methylation states (me1, me2 or me3) 
ptm_ab_pos <- ptm_ab[grepl("K27me", ptm_ab$mi) & grepl("K36me", ptm_ab$mj),]

## Create a crosstalkmap from the filtered data frame that will be saved as a pdf in the 'plots/' folder
CrossTalkMap(
  ptm_ab_pos,
  splitplot_by = "timepoint", # Make one crosstalk map for each time point (Irrelevant in the case of the mESCs dataset that does not contain multiple time points) 
  connected = "tissue", # Visually connect data point accross the value condition 'tissue'
  group_by = "repl", # Group replicates 
  shapecode = "tissue",
  connect_dots = TRUE, # Visually connect data points of a PTM pair accross conditions  
  with_arrows = FALSE, # Use arrow in connection lines  
  colcode = "pij", # Variable to be color encoded ('pi' 'pj' or 'pij')
  which_label = "mimj", # Data points' labels to be diplayed in the crosstalk map ('mj', 'mij')
  contour_lines = TRUE, # Display Interplay score lines 
  contour_labels = "short", 
  col_scheme = "legacy",  
  filename_string = "K27meK36me", 
  filename_ext = "pdf", # Export as pdf
  outdir = "plots/" # output directory 
)





###=- TASK 3.4 -=###
#-Crosstalk maps and PTM abundance/interplay lineplot-#


# Define a list a pairwise combinations of PTMs that will be used to filter the dataframe of PTM abundances.
# For each pair of PTMs a subset dataframe will be created by filtering out rows based on values in the column 'mi'
# matching the first term in a pair and values in the column 'mj' matching the second term of a pair. 
pos_combs <- list(c("K27me1", "K"), c("K27me2", "K"), c("K27me3", "K"))

## TODO explain what parameters connected, group_by and splitplot_by mean here
## (this is another view on the data)
par_connected <- "tissue"
par_group_by <- "repl"
par_splitplot_by <- "timepoint"

# multiple mimj pairs per plot
for (pos_comb in pos_combs) {
  ptm_ab_pos_mod <-
    ptm_ab[grepl(pos_comb[[1]], ptm_ab$mi) &
             grepl(pos_comb[[2]], ptm_ab$mj), ]
  
  
  
  # Show all
  CrossTalkMap(
    ptm_ab_pos_mod,
    splitplot_by = par_splitplot_by,
    connected = par_connected,
    group_by = par_group_by,
    connect_dots = TRUE,
    with_arrows = FALSE,
    colcode = "pj",
    shapecode = "tissue",
    which_label = "mj",
    contour_lines = TRUE,
    contour_labels = "short",
    col_scheme = "standard",
    filename_string = paste0(paste0(pos_comb, collapse = '-'), "_", "all"),
    outdir = "plots/"
  )
  
  
  # Find 10 most abundant m_j and subset data frame
  ptm_ab_mi_ordered <- ptm_ab_pos_mod[order(-ptm_ab_pos_mod$pj), ]
  ptm_ab_mi_uniqmj <- ptm_ab_mi_ordered[!duplicated(ptm_ab_mi_ordered$mj), ]
  ptm_ab_mi_topab <- head(ptm_ab_mi_uniqmj, n = 10)
  topab_mj <- ptm_ab_mi_topab[, "mj"]
  ptm_ab_pos_mod_top <- ptm_ab_pos_mod[grep(paste(topab_mj, collapse = "|"), ptm_ab_pos_mod$mj),]
  
  # Retain data for WT and Suz12-/- conditions only 
  ptm_ab_pos_mod_top_cond <- ptm_ab_pos_mod_top[grepl(c("mESCs|mESCs Suz12-/-"), ptm_ab_pos_mod_top$tissue), ]
  
  # Plot top 10 most abundant PTM pair
  CrossTalkMap(
    ptm_ab_pos_mod_top_cond,
    splitplot_by = par_splitplot_by,
    connected = par_connected,
    group_by = par_group_by,
    connect_dots = TRUE,
    with_arrows = FALSE,
    colcode = "pj",
    shapecode = "tissue",
    which_label = "mj",
    contour_lines = TRUE,
    contour_labels = "short",
    col_scheme = "standard",
    filename_string = paste0(paste0(pos_comb, collapse = '-'), "_", "top10"),
    outdir = "plots/"
  )
  
  
  ## line plots for abundances, co-occurrence, interplay score for each tissue, each top 10 mimj combination
  mis <- unique(ptm_ab_pos_mod$mi)
  
  for (mi in mis) {
    ptm_ab_mi <- ptm_ab_pos_mod[ptm_ab_pos_mod$mi == mi, ]
    
    for (mj in unique(ptm_ab_mi$mj)) {
      ptm_ab_mimj <- ptm_ab_mi[ptm_ab_mi$mj == mj, ]
      
      for (feature in unique(ptm_ab_mimj$timepoint)) {
        ptm_ab_mimj_tis <- ptm_ab_mimj[ptm_ab_mimj$timepoint == feature, ]
        line_ct(
          ptm_ab_mimj_tis,
          connected = par_connected,
          label = "knockouts",
          outdir = "plots/"
        )
        
      }
    }
  }
}
