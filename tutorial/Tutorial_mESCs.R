###########################################################
##             PTM-CrossTalkMapperTutorial               ##
## Analysis of mouse embrionyc stem cells histone 3 PTMs ##
###########################################################


###=- SECTION 1.1 -=###
#-Packages installation and loading-#

##Vector of required packages
packages <- c("broom", "tidyr", "gplots", "gtools", "gtools", "ggplot2", "scales", "metR", "forcats", "ggrepel", "gridExtra", "dplyr")

##If packages are missing run following lines
#install.packages(packages, dependencies = TRUE)
#update.packages(ask = FALSE)

#Load Packages
invisible(lapply(packages, library, character.only = TRUE))



###=- SECTION 1.2 -=###
#-Set working directory and load crosstalk mapper-#

#Define the working directory, by default is set to the location of this script:
setwd(dirname(rstudioapi::getSourceEditorContext()$path)) 
#If it doesn't exist, makes a folder 'plots' in the working directory in which the generated plots will be stored,(note that a warning is given if the folder already exists) 
dir.create(file.path(getwd(), "plots"))
#load PTM-CrosstalkMapper function (If "ptm-crosstalkmapper.R" is not located in the working directory please specify the path to the script)
source("ptm-crosstalkmapper_new.R")



###=- SECTION 2.1 -=###
#-Import dataset and prepare data-#

## TODO how to add other files or get them from crosstalkDB / column organisation and name requirements ? 
mytable <- read.csv("MouseStemCell_CrDB.csv") ###import dataset downloaded from crosstalkDB
colnames(mytable)

## TODO Tell column names need to be adapted
mytable$timepoint <- "None"
mytable$biological.replicate <- mytable$repl <- 1
mytable$tissue <- mytable$cell.type...tissue
write.csv(mytable, "t.csv")

## TODO histvars=TRUE is you wish to distinguish between H3.1 and H3.3. Else set to FALSE
data <- prepPTMdata(csv = "t.csv", histvars = TRUE, avrepls = FALSE)



###=- SECTION 2.1 -=###
#-Compute PTM abundance and interplay score-#

ptm_ab <- calcPTMab(data)
ptm_ab$tissue <- factor(ptm_ab$tissue, levels = c("mESCs","ESCs - Dnmt TKO", "ESCs - Ring1A1B-/-", "mESCs Suz12-/-"))



###=- SECTION 3.1 -=###
#-Bar plot of single modification abundance-#

## TODO explain plot type
#### Number of PTM quantifications
# This function calculate single PTM abundances
plotdat <- unique(ptm_ab[,c("hist","tissue","timepoint","repl","mi","pi")])
# list of PTMs with abundance above threshold (0.1 in this case)
filtered_ptms <- names(which(by(ptm_ab$pi, ptm_ab$mi, max) > 0.1))
# Filtered dataframe
plotdat2 <- plotdat[plotdat$mi %in% filtered_ptms,]
p <- ggplot(plotdat2, aes(fill=tissue, x=reorder(mi, -pi), y=pi)) + 
  geom_bar(position= position_dodge(preserve = "single"), stat="identity") +
  theme_minimal() + 
  theme(text = element_text(size = 15), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  coord_cartesian(clip = "off") + 
  facet_wrap( ~hist) +
  xlab("Modification") 
ggsave(filename = paste0("plots/ptm_distribution.pdf"), plot = p)



###=- SECTION 3.2 -=###
#-Heatmap of single PTMs abundance-#

ptm_ab_hist <- ptm_ab[grepl("H3.1", ptm_ab$hist),]#ptm abundance of one histone variant 

## TODO Explain  idea of heatmaps and what to take from there (single versus double PTMs as information base)
plotdat <- unique(ptm_ab[,c("hist","tissue","timepoint","repl","mi","pi")])
flat_matrix <- reshape2::dcast(plotdat, formula= mi ~ tissue + hist + timepoint + repl, value.var = "pi")

heatmap_all(flat_matrix, showSidebar = "tissue", hscale = "none", 
            title_of_plot = "H3 variants, single PTM frequency", label="single_ptms", outdir = "plots/")


#-heatmap of Co-occurring modification frequencies-#

plotdat <- unique(ptm_ab_hist[,c("hist","tissue","timepoint","repl","mi","mj","pi","pij","I")])
plotdat$mij <- paste(plotdat$mi,plotdat$mj,sep="")
flat_matrix <- reshape2::dcast(plotdat, formula= mij ~ tissue + hist + timepoint + repl, value.var = "pij")

heatmap_all(flat_matrix, showSidebar = "tissue", hscale = "none", 
            title_of_plot = "H3 variants, double PTM frequency", label="double_ptms", outdir = "plots/")


#-Heatmap of Interplay scores-#

plotdat <- unique(ptm_ab_hist[,c("hist","tissue","timepoint","repl","mi","mj","pi","pij","I")])
plotdat$mij <- paste(plotdat$mi,plotdat$mj,sep="")
flat_matrix <- reshape2::dcast(plotdat, formula= mij ~ tissue + hist + timepoint + repl, value.var = "I")
flat_matrix = flat_matrix[!is.na(flat_matrix[,1]),]

heatmap_all(flat_matrix, showSidebar = "tissue", hscale = "none", 
            title_of_plot = "H3 variants, interplay scores", label = "interplay_scores", outdir = "plots/")


#-Heatmap of peptides abundances#

plotdat <- unique(mytable[,c("tissue","modifications","timepoint","repl","quantification")])
flat_matrix <- reshape2::dcast(plotdat, formula= modifications ~ tissue + timepoint + repl, value.var = "quantification", fun.aggregate = sum)

heatmap_all(flat_matrix, showSidebar = "tissue", hscale = "row", 
            title_of_plot = "H3 peptides", label="peptides", outdir = "plots/",autosize=TRUE)



###=- SECTION 3.3 -=###
#-Crosstalk maps of K27 and K36 methylations-#

ptm_ab_pos <- ptm_ab[ grepl("K27me", ptm_ab$mi) & grepl("K36me", ptm_ab$mj), ]


CrossTalkMap(ptm_ab_pos,
             splitplot_by = "timepoint", connected = "tissue", group_by = "repl",
             connect_dots = TRUE, with_arrows = FALSE, colcode = "pj", which_label = "mimj",
             contour_lines = TRUE, contour_labels = "short", col_scheme="legacy",
             filename_string = "K27meK36me", filename_ext = "pdf",
             outdir = "plots/", shapecode = "tissue")





###=- SECTION 3.4 -=###
#-Crosstalk maps and PTM abundance/interplay lineplot -#

##################################################
## ALL INTERACTIONS BETWEEN 2 SELECTED POSITIONS ##
## BOTH AS WELL AS ONLY ONE SPECIFIC MODIFICATION TYPE ##
#########################################################


## TODO Again one can decide about which combinations are of interest to see respective plots
# selected pair-wise combinations
pos_combs <- list(c("K27me1", "K"), c("K27me2", "K"), c("K27me3", "K"))

## TODO explain what parameters connected, group_by and splitplot_by mean here
## (this is another view on the data)
par_connected <- "tissue"
par_group_by <- "repl"
par_splitplot_by <- "timepoint" 

# multiple mimj pairs per plot
for (pos_comb in pos_combs) {
  
  ptm_ab_pos_mod <- ptm_ab[grepl(pos_comb[[1]], ptm_ab$mi) & grepl(pos_comb[[2]], ptm_ab$mj),]
  #ptm_ab_pos_mod <- ptm_ab_pos_mod[grepl(c("mESCs|mESCs Suz12-/-"), ptm_ab_pos_mod$tissue), ]
  
  
  # Show all
  CrossTalkMap(ptm_ab_pos_mod,
               splitplot_by = par_splitplot_by, connected = par_connected, group_by = par_group_by,
               connect_dots = TRUE, with_arrows = FALSE, colcode = "pj", shapecode = "tissue", 
               which_label = "mj", contour_lines = TRUE, contour_labels = "short", col_scheme="standard",
               filename_string = paste0(paste0(pos_comb, collapse = '-'), "_", "all"), outdir = "plots/")
  
  
  #find 10 most abundant m_j and subset data frame
  ptm_ab_mi_ordered <- ptm_ab_pos_mod[order(-ptm_ab_pos_mod$pj),]
  ptm_ab_mi_uniqmj <- ptm_ab_mi_ordered[!duplicated(ptm_ab_mi_ordered$mj),]
  ptm_ab_mi_topab <- head(ptm_ab_mi_uniqmj, n = 10)
  topab_mj <- ptm_ab_mi_topab[,"mj"]
  ptm_ab_pos_mod_top <- ptm_ab_pos_mod[grep(paste(topab_mj, collapse = "|"), ptm_ab_pos_mod$mj), ]
  
  
  # Plot top 10 most abundant m_j
  CrossTalkMap(ptm_ab_pos_mod_top,
               splitplot_by = par_splitplot_by, connected = par_connected, group_by = par_group_by,
               connect_dots = TRUE, with_arrows = FALSE, colcode = "pj", shapecode = "tissue", 
               which_label = "mj", contour_lines = TRUE, contour_labels = "short", col_scheme="standard",
               filename_string = paste0(paste0(pos_comb, collapse = '-'), "_", "top10"), outdir = "plots/")
  
  
  ## line plots for abundances, co-occurrence, interplay score for each tissue, each top 10 mimj combination
  mis <- unique(ptm_ab_pos_mod$mi)
  
  for (mi in mis) {
    ptm_ab_mi <- ptm_ab_pos_mod[ptm_ab_pos_mod$mi == mi,]
    
    for (mj in unique(ptm_ab_mi$mj)) {
      ptm_ab_mimj <- ptm_ab_mi[ptm_ab_mi$mj == mj,]
      
      for (feature in unique(ptm_ab_mimj$timepoint)) {
        ptm_ab_mimj_tis <- ptm_ab_mimj[ptm_ab_mimj$timepoint == feature,]
        line_ct(ptm_ab_mimj_tis, connected=par_connected, label="knockouts", outdir = "plots/")
        
      }
    }
  }
}




#----------------------------------------------------------------------------------------------------


  

  



################################
## LINE PLOT of ptm abundance ##
################################

## line plots for change of abundance over time for each PTM m_i
mi_dat <- unique(ptm_ab_pos[, c("hist", "tissue", "timepoint", "repl", "mi", "pi")])
for (feature in unique(mi_dat$timepoint)) {
  mi_dat_tis <- mi_dat[mi_dat$timepoint == feature,]
  line_ab(mi_dat_tis, connected="tissue", label="knockouts", outdir = "plots/")
}


## line plots for abundances, co-occurrence, interplay score for each tissue, each mimj combination
for (mj in unique(ptm_ab_pos$mj)) {
  ptm_ab_mimj <- ptm_ab_pos[ptm_ab_pos$mj == mj,]
  for (feature in unique(ptm_ab_mimj$timepoint)) {
    ptm_ab_mimj_tis <- ptm_ab_mimj[ptm_ab_mimj$timepoint == feature,]
    line_ct(ptm_ab_mimj_tis, connected="tissue", label="knockouts", outdir = "plots/")
  }
}



#########################################################
##                   CROSSTALK MAP                     ##
## BOTH AS WELL AS ONLY ONE SPECIFIC MODIFICATION TYPE ##
#########################################################

## ALL / TOP X INTERACTIONS OF 1 SELECTED MODIFICATION ##
## INCL LINE PLOT FOR EACH COMBINATION ##

## one plot per PTM m_i with top x most abundant co-occurring m_j
### Figures for total histones and histone variants

## TODO explain what parameters connected, group_by and splitplot_by mean (is in the documentation)
par_connected <- "tissue"
par_group_by <- "repl"
par_splitplot_by <- "timepoint" 

## TODO select the positions of interest here. Too many will create too many files and too long running time
#### positions of interest ###
poi <- c("K27me")


## TODO explain that this for loop creates various crosstalk maps and line plots that will be placed
## into the "plots" folder
for (mi_pos in poi) {
  
  ptm_ab_pos <- ptm_ab[grep(mi_pos, ptm_ab$mi),]
  mis <- unique(ptm_ab_pos$mi)
  
  for (mi in mis) {
    
    ptm_ab_mi <- ptm_ab_pos[ptm_ab_pos$mi == mi,]
    
    ## no filter
    CrossTalkMap(ptm_ab_mi, splitplot_by = par_splitplot_by, colcode = "pj", connected = par_connected, 
                 group_by = par_group_by, shapecode="tissue",connect_dots = TRUE, with_arrows = FALSE,
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
      CrossTalkMap(midat_plot, splitplot_by = par_splitplot_by, colcode = "pj", connected = par_connected, 
                   group_by = par_group_by, shapecode="tissue",connect_dots = TRUE, with_arrows = FALSE,
                   filename_string = paste0(mi, "_pj-variants-top-", top), outdir = "plots/")
      
      # For total histones irrespective of variants
      midat_plot$hist <- "Total H3"
      CrossTalkMap(midat_plot, splitplot_by = par_splitplot_by, colcode = "pj", connected = par_connected, 
                   group_by = par_group_by, shapecode="tissue",connect_dots = TRUE, with_arrows = FALSE,
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









############function modified################
heatmap_all <- function(flat_matrix, showSidebar = "tissue", hscale="none", title_of_plot = "",label="", outdir = getwd(),autosize=FALSE) {
  ## function to create heatmaps with a sidebar that shows the different categories 
  # showSidebar: category to show with clustered heatmap
  # hscale: used scaling: none, row or column
  # title_of_plot
  # label: label to be added to pdf-file
  # outdir: folder to place pdf-file
  # autosize: automatically modify heatmap layout to display all labels 
  
  rownames(flat_matrix) <- flat_matrix[,1]
  flat_matrix <- flat_matrix[,2:ncol(flat_matrix)]
  flat_matrix <- as.matrix(flat_matrix)
  flat_matrix[!is.finite(flat_matrix)] <- NA
  flat_matrix <- flat_matrix[rowSums(!is.na(flat_matrix)) > ncol(flat_matrix) / 2, ]
  
  # Decide which features will be shown as side bar
  colvec <- vector(,ncol(flat_matrix))
  
  tfeatures <- unique(ptm_ab[[showSidebar]])
  for (i in 1:length(tfeatures)) {
    colvec[grep(tfeatures[i], colnames(flat_matrix))] <- i
  }
  
  #automatically adjust heatmap's layout parameters to display all the row labels
  if(autosize==TRUE){
    pdfHeight = 7+(nrow(flat_matrix)/10) #height of the pdf based on the number of rows in the heatmap
    lheiM = c(2/pdfHeight*10,10) # set the top dendrogram size / heatmap size ratio
    mar = c(15,8)
  }else{ #by default
    pdfHeight = 7
    lheiM = NULL
    mar = c(5,5)
  }
  
  pdf(paste0(outdir, "heatmap_all_", label, ".pdf"), height=pdfHeight)
  
  heatmap.2(as.matrix(flat_matrix), Rowv = T, Colv = T, cexRow = 0.6, cexCol = 1,
            main = title_of_plot, trace = "none",
            scale = hscale, col=bluered(100), density.info = "density", 
            ColSideColors = rainbow(length(colvec))[colvec], srtCol=45, margins=mar, lhei = lheiM)
  
  dev.off()
}


