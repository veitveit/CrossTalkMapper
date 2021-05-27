# assumes working dir to be in mouse-tissue-analysis/
source("../ctm-functions/ptm-crosstalkmapper_old.R")

# Case (1):
# Histone H3 total, averaged replicates

####################
## Specific input ##
####################

# CrosstalkDB data
ctdb_data <- "../data/mouse-tissues_ctdb_timerep_4timepoints.csv"

#########################
## Data Pre-Processing ##
#########################

data <- prepPTMdata(ctdb_data, histvars = FALSE, avrepls = TRUE)

# shorten tissue labels, specific for dataset PXD005300
data$cell.type...tissue <- gsub(".*, ", "", data$cell.type...tissue)
data$cell.type...tissue <- paste0(toupper(substr(data$cell.type...tissue, 1, 1)),
                                  substr(data$cell.type...tissue, 2, nchar(data$cell.type...tissue)))

###############################
## PTM abundance calculation ##
###############################

ptm_ab <- calcPTMab(data, outdir = "h3/data/")

#####################
## Filter and plot ##
#####################

# positions of interest
poi <- c("K9", "K14", "K18", "K23", "K27", "K36")

## ALL / TOP X INTERACTIONS OF 1 SELECTED MODIFICATION ##
## INCL LINE PLOT FOR EACH COMBINATION ##

# one plot per PTM m_i with top x most abundant co-occurring m_j
for (mi_pos in poi) {
  
  ptm_ab_pos <- ptm_ab[grep(mi_pos, ptm_ab$mi),]
  mis <- unique(ptm_ab_pos$mi)
  
  for (mi in mis) {
    
    ptm_ab_mi <- ptm_ab_pos[ptm_ab_pos$mi == mi,]
    
    ## no filter
    CrossTalkMap(ptm_ab_mi, splitplot_by = "tissue", colcode = "pj", connected = "timepoint", group_by = "repl",
                     connect_dots = TRUE, with_arrows = TRUE,
                     filename_string = paste0(mi, "_mj-all"), outdir = "h3/crosstalkmaps/")

    ## top 10/5 abundant m_j
    # rank m_i data by abundance of m_j, regardless of time point, tissue or histone variant
    ptm_ab_mi_ordered <- ptm_ab_mi[order(-ptm_ab_mi$pj),]
    ptm_ab_mi_uniqmj <- ptm_ab_mi_ordered[!duplicated(ptm_ab_mi_ordered$mj),]
    # find 5 and 10 most abundant m_j, subset data frame and plot only these data points
    for (top in c(5, 10)) {
      ptm_ab_mi_topab <- head(ptm_ab_mi_uniqmj, n = top)
      topab_mj <- ptm_ab_mi_topab[,"mj"]

      midat_plot <- ptm_ab_mi[grep(paste(topab_mj, collapse = "|"), ptm_ab_mi$mj),]

      CrossTalkMap(midat_plot, splitplot_by = "tissue", colcode = "pj", connected = "timepoint", group_by = "repl",
                       connect_dots = TRUE, with_arrows = TRUE,
                       filename_string = paste0(mi, "_pj-top-", top), outdir = "h3/crosstalkmaps/")
    }
    
    ## line plots for change of abundance over time for each PTM m_i
    mi_dat <- unique(ptm_ab_mi[, c("hist", "tissue", "timepoint", "repl", "mi", "pi")])
    for (tissue in unique(mi_dat$tissue)) {
      mi_dat_tis <- mi_dat[mi_dat$tissue == tissue,]
      line_ab(mi_dat_tis, outdir = "h3/lineplots_pi/")
    }
    
    ## line plots for abundances, co-occurrence, interplay score for each tissue, each mimj combination
    for (mj in unique(ptm_ab_mi$mj)) {
      ptm_ab_mimj <- ptm_ab_mi[ptm_ab_mi$mj == mj,]
      for (tissue in unique(ptm_ab_mimj$tissue)[1]) {
        ptm_ab_mimj_tis <- ptm_ab_mimj[ptm_ab_mimj$tissue == tissue,]
        line_ct(ptm_ab_mimj_tis, outdir = "h3/lineplots_ct-params/")
      }
    }
    
  }
  
}

## ALL INTERACTIONS BETWEEN 2 SELECTED POSITIONS ##
## BOTH AS WELL AS ONLY ONE SPECIFIC MODIFICATION TYPE ##

# selected pair-wise combinations
pos_combs <- list(c("K9", "K27"), c("K9", "K36"), c("K27", "K36"), c("K14", "K18"), c("K14", "K23"), c("K18", "K23"))

# multiple mimj pairs per plot
# but only 2 positions per (multi) plot
for (pos_comb in pos_combs) {
  ptm_ab_pos <- ptm_ab[grepl(pos_comb[[1]], ptm_ab$mi) & grepl(pos_comb[[2]], ptm_ab$mj),]
  CrossTalkMap(ptm_ab_pos, splitplot_by = "tissue", colcode = "pj", connected = "timepoint", group_by = "repl",
                   connect_dots = TRUE, with_arrows = TRUE, which_label = "mimj",
                   filename_string = paste0(paste0(pos_comb, collapse = '-'), "_all-ptms"), outdir = "h3/crosstalkmaps/")
  # only show methylations
  mod <- "me"
  ptm_ab_pos_mod <- ptm_ab_pos[grepl(mod, ptm_ab_pos$mi) & grepl(mod, ptm_ab_pos$mj),]
  if (nrow(ptm_ab_pos_mod) > 0) {
    CrossTalkMap(ptm_ab_pos_mod, splitplot_by = "tissue", colcode = "pj", connected = "timepoint", group_by = "repl",
                     connect_dots = TRUE, with_arrows = TRUE, which_label = "mimj",
                     filename_string = paste0(paste0(pos_comb, collapse = '-'), "_", mod), outdir = "h3/crosstalkmaps/")
  }
}