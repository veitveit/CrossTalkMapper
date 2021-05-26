# assumes working dir to be in mouse-tissue-analysis/
source("../ctm-functions/ptm-crosstalkmapper.R")

# Case (3):
# Histone H3 total, individual replicates

####################
## Specific input ##
####################

# CrosstalkDB data
ctdb_data <- "../data/mouse-tissues_ctdb_timerep_4timepoints.csv"

# Output dirs
ab_out_dir <- "repls/data/"
ctm_out_dir <- "repls/crosstalkmaps/"
# ab_out_dir <- ctm_out_dir <- "test/"

#########################
## Data Pre-Processing ##
#########################

data <- prepPTMdata(ctdb_data, histvars = FALSE, avrepls = FALSE)

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
poi <- c("K9", "K27", "K36")

## TOP X INTERACTIONS OF 1 SELECTED MODIFICATION ##

# one plot per PTM m_i with top x most abundant co-occurring m_j, individual replicates
for (mi_pos in poi[1]) {
  
  ptm_ab_pos <- ptm_ab[grep(mi_pos, ptm_ab$mi),]
  mis <- unique(ptm_ab_pos$mi)
  
  for (mi in mis[1]) {

    ptm_ab_mi <- ptm_ab_pos[ptm_ab_pos$mi == mi,]
    
    ## top 2/3/5 abundant m_j
    # rank m_i data by abundance of m_j, regardless of time point, tissue or histone variant
    ptm_ab_mi_ordered <- ptm_ab_mi[order(-ptm_ab_mi$pj),]
    ptm_ab_mi_uniqmj <- ptm_ab_mi_ordered[!duplicated(ptm_ab_mi_ordered$mj),]
    # find 2, 3 and 5 most abundant m_j, subset data frame and plot only these data points
    for (top in c(2, 3, 5)[1]) {
      ptm_ab_mi_topab <- head(ptm_ab_mi_uniqmj, n = top)
      topab_mj <- ptm_ab_mi_topab[,"mj"]
      
      midat_plot <- ptm_ab_mi[grep(paste(topab_mj, collapse = "|"), ptm_ab_mi$mj),]
      
      CrossTalkMap(midat_plot, splitplot_by = "timepoint", colcode = "tissue", connected = "repl", group_by = "tissue",
                       connect_dots = TRUE, with_arrows = FALSE, hide_axes = TRUE)#,
                       #filename_string = paste0(mi, "_pj-top-", top), outdir = ctm_out_dir)
        
      }

  #   ## plot only combinations of PTMs at positions of interest
  #   ptm_ab_poi <- ptm_ab_mi[grep(paste0(poi, collapse = "|"), ptm_ab_mi$mj),]
  #   CrossTalkMap(ptm_ab_poi, filename_string = paste0(poi, collapse = ""), outdir = "plots/", indplots = FALSE)
  #     
  }
  
}
