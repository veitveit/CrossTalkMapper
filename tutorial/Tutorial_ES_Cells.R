## TODO comment about installation
# install.packages(c("broom", "tidyr", "gplots", "gtools", "gtools", "ggplot2", "scales", "metR", "forcats", "ggrepel", "gridExtra"), dependencies = TRUE)
# update.packages(ask = FALSE)

library(broom)
library(gplots)
library(tidyr)
library(tools)
library(gtools)
library(ggplot2)
library(scales)
library(metR)
library(ggrepel)
library(gridExtra)
library(forcats)
library(stringr)
library(reshape2)

## TODO set folder
setwd("./")
source("ptm-crosstalkmapper_new.R")

## TODO how to add other files or get them from crosstalkdb
mytable <- read.csv("MouseStemCells.csv") #check the data has correct format
colnames(mytable)

## TODO Tell column names need to be adappted
mytable$timepoint <- "None"
mytable$biological.replicate <- mytable$repl <- 1
mytable$tissue <- mytable$cell.type...tissue
write.csv(mytable, "t.csv")


## TODO histvars=TRUE is you wish to distinguish between H3.1 and H3.3. Else set to FALSE
data <- prepPTMdata(csv = "t.csv", histvars = TRUE, avrepls = FALSE)
ptm_ab <- calcPTMab(data)
head(ptm_ab)

## TODO: explain how to change example
# Example for K9 and K27 PTMs
ptm_ab_pos <- ptm_ab[grepl("K9", ptm_ab$mi) & grepl("K27", ptm_ab$mj),]
# TODO create folder "plots" in working directory
CrossTalkMap(ptm_ab_pos,
             splitplot_by = "timepoint", connected = "tissue", group_by = "repl",
             connect_dots = TRUE, with_arrows = FALSE, colcode = "pj", which_label = "mimj",
             contour_lines = TRUE, contour_labels = "short",
             filename_string = "K9-K27_all-ptms", filename_ext = "pdf",
             outdir = "plots/")

## TODO select the positions of interest here. Too many will create too many files and too long running time
#### positions of interest ###
poi <- c("K9", "K27", "K36")

## ALL / TOP X INTERACTIONS OF 1 SELECTED MODIFICATION ##
## INCL LINE PLOT FOR EACH COMBINATION ##

## one plot per PTM m_i with top x most abundant co-occurring m_j
### Figures for total histones and histone variants

## TODO explain what parameters connected, group_by and splitplot_by mean (is in the documentation)
par_connected <- "tissue"
par_group_by <- "repl"
par_splitplot_by <- "timepoint" 


## TODO explain that this for loop creates various crosstalk maps and line plots that will be placed
## into the "plots" folder
for (mi_pos in poi) {
  
  ptm_ab_pos <- ptm_ab[grep(mi_pos, ptm_ab$mi),]
  mis <- unique(ptm_ab_pos$mi)
  
  for (mi in mis) {
    
    ptm_ab_mi <- ptm_ab_pos[ptm_ab_pos$mi == mi,]
    
    ## no filter
    CrossTalkMap(ptm_ab_mi, splitplot_by = par_splitplot_by, colcode = "pj", connected = par_connected, group_by = par_group_by,
                 connect_dots = TRUE, with_arrows = TRUE,
                 filename_string = paste0(mi, "_mj-all"), outdir = "plots/")
    
    ## top 10/5 abundant m_j
    # rank m_i data by abundance of m_j, regardless of time point, tissue or histone variant
    ptm_ab_mi_ordered <- ptm_ab_mi[order(-ptm_ab_mi$pj),]
    ptm_ab_mi_uniqmj <- ptm_ab_mi_ordered[!duplicated(ptm_ab_mi_ordered$mj),]
    # find 5 and 10 most abundant m_j, subset data frame and plot only these data points
    for (top in c(5, 10)) {
      ptm_ab_mi_topab <- head(ptm_ab_mi_uniqmj, n = top)
      topab_mj <- ptm_ab_mi_topab[,"mj"]
      
      midat_plot <- ptm_ab_mi[grep(paste(topab_mj, collapse = "|"), ptm_ab_mi$mj),]
      
      # Separate figures for variants
      CrossTalkMap(midat_plot, splitplot_by = par_splitplot_by, colcode = "pj", connected = par_connected, group_by = par_group_by,
                   connect_dots = TRUE, with_arrows = TRUE,
                   filename_string = paste0(mi, "_pj-variants-top-", top), outdir = "plots/")
      # For total histones irrespective of variants
      midat_plot$hist <- "Total H3"
      CrossTalkMap(midat_plot, splitplot_by = par_splitplot_by, colcode = "pj", connected = par_connected, group_by = par_group_by,
                   connect_dots = TRUE, with_arrows = TRUE,
                   filename_string = paste0(mi, "_pj-total-top-", top), outdir = "plots/")
      
    }
    
    ## line plots for change of abundance over time for each PTM m_i
    mi_dat <- unique(ptm_ab_mi[, c("hist", "tissue", "timepoint", "repl", "mi", "pi")])
    for (feature in unique(mi_dat$timepoint)) {
      mi_dat_tis <- mi_dat[mi_dat$timepoint == feature,]
      line_ab(mi_dat_tis, connected=par_connected, label="knockouts", outdir = "plots/")
    }
    
    ## line plots for abundances, co-occurrence, interplay score for each tissue, each mimj combination
    for (mj in unique(ptm_ab_mi$mj)) {
      ptm_ab_mimj <- ptm_ab_mi[ptm_ab_mi$mj == mj,]
      for (feature in unique(ptm_ab_mimj$timepoint)) {
        ptm_ab_mimj_tis <- ptm_ab_mimj[ptm_ab_mimj$timepoint == feature,]
        line_ct(ptm_ab_mimj_tis, connected=par_connected, label="knockouts", outdir = "plots/")
      }
    }
    
  }
  
}


## ALL INTERACTIONS BETWEEN 2 SELECTED POSITIONS ##
## BOTH AS WELL AS ONLY ONE SPECIFIC MODIFICATION TYPE ##

## TODO Again one can decide about which combinations are of interest to see respective plots
# selected pair-wise combinations
pos_combs <- list(c("K9", "K27"), c("K9", "K36"), 
                  c("K14", "K18"), c("K14", "K27"), c("K14", "K36"),
                  c("K18", "K23"), c("K18", "K27"), c("K18", "K36"), 
                  c("K23", "K27"), c("K23", "K36"), c("K27", "K36"))

## TODO explain what parameters connected, group_by and splitplot_by mean here
## (this is another view on the data)
par_connected <- "timepoint"
par_group_by <- "repl"
par_splitplot_by <- "tissue" 

# multiple mimj pairs per plot
# but only 2 positions per (multi) plot
for (pos_comb in pos_combs) {
  ptm_ab_pos <- ptm_ab[grepl(pos_comb[[1]], ptm_ab$mi) & grepl(pos_comb[[2]], ptm_ab$mj),]
  ptm_ab_pos_mod <- rbind(ptm_ab_pos[grepl("ac", ptm_ab_pos$mi) & grepl("ac", ptm_ab_pos$mj),],
                          ptm_ab_pos[grepl("ac", ptm_ab_pos$mi) & grepl("me", ptm_ab_pos$mj),],
                          ptm_ab_pos[grepl("me", ptm_ab_pos$mi) & grepl("ac", ptm_ab_pos$mj),])
  
  # only show methylations
  mod <- "all"
  CrossTalkMap(ptm_ab_pos_mod, splitplot_by = par_splitplot_by, colcode = "pj", connected = par_connected, group_by = par_group_by,
               connect_dots = TRUE, with_arrows = TRUE, which_label = "mimj",
               filename_string = paste0(paste0(pos_comb, collapse = '-'), "_", mod), outdir = "plots/")
  # only show methylations
  mod <- "me"
  ptm_ab_pos_mod <- ptm_ab_pos[grepl(mod, ptm_ab_pos$mi) & grepl(mod, ptm_ab_pos$mj),]
  if (nrow(ptm_ab_pos_mod) > 0) {
    CrossTalkMap(ptm_ab_pos_mod, splitplot_by = par_splitplot_by, colcode = "pj", connected = par_connected, group_by = par_group_by,
                 connect_dots = TRUE, with_arrows = TRUE, which_label = "mimj",
                 filename_string = paste0(paste0(pos_comb, collapse = '-'), "_", mod), outdir = "plots/")
  }
}


##################################
##CREATE SINGLE MOD BAR PLOTS H4##
##################################

## TODO explain plot type
#### Number of PTM quantifications
## Function to calculate single PTM abundances
plotdat <- unique(ptm_ab[,c("hist","tissue","timepoint","repl","mi","pi")])
# everything with highest value above threshold
filtered_ptms <- names(which(by(ptm_ab$pi, ptm_ab$mi, max) > 0.1))
plotdat2 <- plotdat[plotdat$mi %in% filtered_ptms,]
p <- ggplot(plotdat2, aes(fill=tissue, x=mi, y=pi)) + 
  geom_bar(position= position_dodge(preserve = "single"), stat="identity") +
  theme_minimal() + 
  theme(text = element_text(size = 15), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  coord_cartesian(clip = "off") + 
  facet_wrap( ~hist) 
ggsave(filename = paste0("plots/ptm_distribution.pdf"), plot = p)


###############################################
#### heatmaps of single PTMs
###############################################


## TODO Explain idea of heatmaps and what to take from there (single versus double PTMs as information base)
plotdat <- unique(ptm_ab[,c("hist","tissue","timepoint","repl","mi","pi")])
flat_matrix <- reshape2::dcast(plotdat, formula= mi ~ tissue + hist + timepoint + repl, value.var = "pi")
heatmap_all(flat_matrix, showSidebar = "tissue", hscale = "none", 
            title_of_plot = "H3 variants, single PTM frequency", label="single_ptms", outdir = "plots/")


####################################################
## HEATMAP OF COOCCURING MODIFICATION FREQUENCIES ##
####################################################

plotdat <- unique(ptm_ab[,c("hist","tissue","timepoint","repl","mi","mj","pi","pij","I")])
plotdat$mij <- paste(plotdat$mi,plotdat$mj,sep="")

flat_matrix <- reshape2::dcast(plotdat, formula= mij ~ tissue + hist + timepoint + repl, value.var = "pij")
heatmap_all(flat_matrix, showSidebar = "tissue", hscale = "none", 
            title_of_plot = "H3 variants, double PTM frequency", label="double_ptms", outdir = "plots/")


####################################################
## HEATMAP OF Interplay scores ##
####################################################


plotdat <- unique(ptm_ab[,c("hist","tissue","timepoint","repl","mi","mj","pi","pij","I")])
plotdat$mij <- paste(plotdat$mi,plotdat$mj,sep="")

flat_matrix <- reshape2::dcast(plotdat, formula= mij ~ tissue + hist + timepoint + repl, value.var = "I")
heatmap_all(flat_matrix, showSidebar = "tissue", hscale = "none", 
            title_of_plot = "H3 variants, interplay scores", label="interplay_scores", outdir = "plots/")



####################################################
## HEATMAP OF peptides ##
####################################################


plotdat <- unique(mytable[,c("tissue","modifications","timepoint","repl","quantification")])
flat_matrix <- reshape2::dcast(plotdat, formula= modifications ~ tissue + timepoint + repl, value.var = "quantification", fun.aggregate = sum)
heatmap_all(flat_matrix, showSidebar = "tissue", hscale = "row", 
            title_of_plot = "H3 peptides", label="peptides", outdir = "plots/")


