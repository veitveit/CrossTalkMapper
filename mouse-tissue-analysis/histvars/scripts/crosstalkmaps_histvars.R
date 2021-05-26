# assumes working dir to be in mouse-tissue-analysis/
source("../ctm-functions/ptm-crosstalkmapper.R")

# Case (2):
# Histone variants treated individually, averaged replicates

####################
## Specific input ##
####################

# CrosstalkDB data
ctdb_data <- "../data/mouse-tissues_ctdb_timerep_4timepoints.csv"

# Output dirs
ab_out_dir <- "histvars/data/"
ctm_out_dir <- "histvars/crosstalkmaps/"
line_ab_out_dir <- "histvars/lineplots_pi/"
line_ct_out_dir <- "histvars/lineplots_ct-params/"
# ab_out_dir <- ctm_out_dir <- line_ab_out_dir <- line_ct_out_dir <- "test/"

#########################
## Data Pre-Processing ##
#########################

data <- prepPTMdata(ctdb_data, histvars = TRUE, avrepls = TRUE)

# shorten tissue labels, specific for dataset PXD005300
data$cell.type...tissue <- gsub(".*, ", "", data$cell.type...tissue)
data$cell.type...tissue <- paste0(toupper(substr(data$cell.type...tissue, 1, 1)),
                                  substr(data$cell.type...tissue, 2, nchar(data$cell.type...tissue)))

###############################
## PTM abundance calculation ##
###############################

ptm_ab <- calcPTMab(data, outdir = ab_out_dir)

#####################
## Filter and plot ##
#####################

# positions of interest
poi <- c("K9", "K14", "K18", "K23", "R26", "K27", "K36")

## ALL / TOP X INTERACTIONS OF 1 SELECTED MODIFICATION ##
## INCL LINE PLOT FOR EACH COMBINATION ##

# one plot per PTM m_i with top x most abundant co-occurring m_j
for (mi_pos in poi) {
  
  ptm_ab_pos <- ptm_ab[grep(mi_pos, ptm_ab$mi),]
  mis <- unique(ptm_ab_pos$mi)
  
  for (mi in mis) {
    
    ptm_ab_mi <- ptm_ab_pos[ptm_ab_pos$mi == mi,]
    
    # no filter
    CrossTalkMap(ptm_ab_mi, splitplot_by = "tissue", colcode = "pj", connected = "timepoint", group_by = "repl",
                     connect_dots = TRUE, with_arrows = TRUE,
                     filename_string = paste0(mi, "_mj-all"), outdir = ctm_out_dir)

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
                       filename_string = paste0(mi, "_pj-top-", top), outdir = ctm_out_dir)
    }
    
    ## line plots for change of abundance over time for each PTM m_i and histone variant
    mi_dat <- unique(ptm_ab_mi[, c("hist", "tissue", "timepoint", "repl", "mi", "pi")])
    for (tissue in unique(mi_dat$tissue)) {
      mi_dat_tis <- mi_dat[mi_dat$tissue == tissue,]
      for (hist in unique(mi_dat_tis$hist)) {
        mi_dat_tis_hist <- mi_dat_tis[mi_dat_tis$hist == hist,]
        line_ab(mi_dat_tis_hist, outdir = line_ab_out_dir)
      }
    }
    
    ## line plots for abundances, co-occurrence, interplay score for each tissue, histone variant and mimj combination
    for (mj in unique(ptm_ab_mi$mj)) {
      ptm_ab_mimj <- ptm_ab_mi[ptm_ab_mi$mj == mj,]
      for (tissue in unique(ptm_ab_mimj$tissue)) {
        ptm_ab_mimj_tis <- ptm_ab_mimj[ptm_ab_mimj$tissue == tissue,]
        for (hist in unique(ptm_ab_mimj_tis$hist)) {
          ptm_ab_mimj_tis_hist <- ptm_ab_mimj_tis[ptm_ab_mimj_tis$hist == hist,]
          line_ct(ptm_ab_mimj_tis_hist, outdir = line_ct_out_dir)
        }
      }
    }

  }
  
}

## ALL INTERACTIONS BETWEEN 2 SELECTED POSITIONS ##
## BOTH AS WELL AS ONLY ONE SPECIFIC MODIFICATION TYPE ##

## selected pair-wise combinations, one per plot ##

# selected pair-wise combinations
pos_combs <- list(c("K9", "K27"), c("K9", "K36"), c("K27", "K36"), c("K14", "K18"), c("K14", "K23"), c("K18", "K23"))

# multiple mimj pairs per plot
# but only 2 positions per (multi) plot
for (pos_comb in pos_combs) {
  ptm_ab_pos <- ptm_ab[grepl(pos_comb[[1]], ptm_ab$mi) & grepl(pos_comb[[2]], ptm_ab$mj),]
  CrossTalkMap(ptm_ab_pos, splitplot_by = "tissue", colcode = "pj", connected = "timepoint", group_by = "repl",
                   connect_dots = TRUE, with_arrows = TRUE, which_label = "mimj",
                   filename_string = paste0(paste0(pos_comb, collapse = '-'), "_all-ptms"), outdir = ctm_out_dir)
  # only show methylations
  mod <- "me"
  ptm_ab_pos_mod <- ptm_ab_pos[grepl(mod, ptm_ab_pos$mi) & grepl(mod, ptm_ab_pos$mj),]
  if (nrow(ptm_ab_pos_mod) > 0) {
    CrossTalkMap(ptm_ab_pos_mod, splitplot_by = "tissue", colcode = "pj", connected = "timepoint", group_by = "repl",
                     connect_dots = TRUE, with_arrows = TRUE, which_label = "mimj",
                     filename_string = paste0(paste0(pos_comb, collapse = '-'), "_", mod), outdir = ctm_out_dir)
  }
}

## ALL INTERACTIONS BETWEEN ALL SELECTED POSITIONS ##
## ONLY FOR ACETYLATION ##

## all pair-wise ac combinations in one plot ##

# generate all pair-wise combinations of positions
# if positions in list are ordered numerically (K9 < R26 < K27 < K36 etc),
# the PTMs within the combinations are in the right order as well
pos_combs <- combn(poi, 2)
mod <- "ac"

# multiple mimj pairs per plot, all positions in one plot (ac only)
# extract all relevant data entries
ptm_ab_pos <- data.frame()
for (pos_comb in seq(1, ncol(pos_combs))) {
  ptm_ab_pos <- rbind(ptm_ab_pos,
                      ptm_ab[ptm_ab$mi == paste0(pos_combs[1,pos_comb], mod) &
                               ptm_ab$mj == paste0(pos_combs[2,pos_comb], mod),])
}
CrossTalkMap(ptm_ab_pos, splitplot_by = "tissue", colcode = "pj", connected = "timepoint", group_by = "repl",
                 connect_dots = TRUE, with_arrows = TRUE, which_label = "mimj",
                 filename_string = paste0(paste0(poi, collapse = '-'), "_all-ac"), outdir = ctm_out_dir)
