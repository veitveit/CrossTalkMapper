# assumes working dir to be in mouse-tissue-analysis/
source("../ctm-functions/ptm-crosstalkmapper.R")

# Case (1):
# Histone H3 total, averaged replicates

####################
## Specific input ##
####################

# CrosstalkDB data
ctdb_data <- "../data/mouse-tissues_ctdb_timerep_4timepoints.csv"

# Output dirs
ab_out_dir <- "h3/data/"
ctm_out_dir <- "h3/crosstalkmaps/"
line_ab_out_dir <- "h3/lineplots_pi/"
line_ct_out_dir <- "h3/lineplots_ct-params/"
# ab_out_dir <- ctm_out_dir <- line_ab_out_dir <- line_ct_out_dir <- "test/"

#######################################
## Data Retrieval and Pre-Processing ##
#######################################

## from CrosstalkDB - two possibilities

# define range of CrosstalkDB identifier integers (e.g. 13 for CrDB000013)
crdb_ints <- seq(62,129)
## filtering should not be necessary when non-existing IDs are handled in download
# # define integers that should be removed from that range (do not exist) - not necessary, then download is just skipped
# rm_ints <- c(seq(64,67))
# # only keep the remaining integers
# crdb_ints <- crdb_ints_range[! crdb_ints_range %in% rm_ints]
# convert integers to CrDB identifiers
# crdb_ids <- gen_CrDB_ids(crdb_ints)
# convert integers to CrDB identifiers --> done in prep_CrDBdata()
# crdb_ids <- gen_CrDB_ids(crdb_ints_range)
# retrieve datasets from CrosstalkDB
crdb_data <- prep_CrDBdata(crdb_ints,
                           out_file = "~/Arbeit/prgroup-sdu/projects/CrossTalkMapper-dev/test/crdb_data.tsv",
                           ## DEBUG
                           dl_test = "new_descr")

# filter dataset (exclude 24 months samples)
# can also be done before on identifier level, but is easier here by excluding everything "24 months"
# (no need to know the specific identifiers of 24-mo samples)


# alternatively, supply a file of CrDB identifiers and conditions to prep_CrDBdata()
cond_file <- "~/Arbeit/prgroup-sdu/projects/CrossTalkMapper-dev/data/sample_conditions_no24mo.tsv"
crdb_data <- prep_CrDBdata(cond_file = cond_file, dl_test = "old")


## alternatively, get PTM data from file with CrDB data file structure (no quotes required)
ctdb_data <- "../data/mouse-tissues_ctdb_timerep_4timepoints.csv"

#########################
## PTM Data Processing ##
#########################

data <- prepPTMdata(ctdb_data, histvars = FALSE, avrepls = TRUE)

# # shouldn't be necessary anymore
# # shorten tissue labels, specific for dataset PXD005300
# data$cell.type...tissue <- gsub(".*, ", "", data$cell.type...tissue)
# data$cell.type...tissue <- paste0(toupper(substr(data$cell.type...tissue, 1, 1)),
#                                   substr(data$cell.type...tissue, 2, nchar(data$cell.type...tissue)))

###############################
## PTM abundance calculation ##
###############################

ptm_ab <- calcPTMab(data, outdir = ab_out_dir)

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
    
    ## line plots for change of abundance over time for each PTM m_i
    mi_dat <- unique(ptm_ab_mi[, c("hist", "tissue", "timepoint", "repl", "mi", "pi")])
    for (tissue in unique(mi_dat$tissue)) {
      mi_dat_tis <- mi_dat[mi_dat$tissue == tissue,]
      line_ab(mi_dat_tis, outdir = line_ab_out_dir)
    }
    
    ## line plots for abundances, co-occurrence, interplay score for each tissue, each mimj combination
    for (mj in unique(ptm_ab_mi$mj)) {
      ptm_ab_mimj <- ptm_ab_mi[ptm_ab_mi$mj == mj,]
      for (tissue in unique(ptm_ab_mimj$tissue)[1]) {
        ptm_ab_mimj_tis <- ptm_ab_mimj[ptm_ab_mimj$tissue == tissue,]
        line_ct(ptm_ab_mimj_tis, outdir = line_ct_out_dir)
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