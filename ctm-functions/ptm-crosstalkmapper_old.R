library(tidyr)
library(tools)
library(ggplot2)
library(scales)
library(metR)
library(ggrepel)
library(gridExtra)

# sample = one combination of tissue + timepoint/any conditions (+ replicate before averaging)
# for each sample, all quantifications have to add up to 1

#########################
## Data Pre-processing ##
#########################

cleanup_ctdb <- function(data) {
  # remove unnecessary columns from CrosstalkDB-derived .csv files
  data[, c("X", "accession.number", "species", "number.of.modifications", "date", "user", "MS.resolution", "Link.to.original.data",
             "peptide.Sequence", "dataset.id")] <- list(NULL)
  return(data)
}

zero_imp <- function(data) {
  # 0 imputation
  # afterwards, every sample-protein combination has the same set of modifications and proteins
  # i.e. for each modification contained in the dataset there is one observation / row per sample-protein combination
  # (quantification = 0 if not measured)
  # only combinations of tissue, timepoint and replicate that actually exist in the dataset are included
  # (relevant in case of differing samples numbers for different conditions, e.g. fewer time points for a subset of tissues)
  data_0 <- complete(data, nesting(cell.type...tissue, timepoint, biological.replicate), protein.name, modifications,
                     fill = list(quantification = 0))
  return(data_0)
}

av_repls <- function(data) {
  # average replicates
  avdata <- aggregate(data$quantification, by = data[c("cell.type...tissue", "timepoint", "protein.name", "modifications")], mean)
  names(avdata)[names(avdata) == 'x'] <- "quantification"
  avdata$biological.replicate <- "av"
  return(avdata)
}

histvarquant <- function(data) {
  # calculate histone variant abundances for each combination of conditions
  histvarquants <- data.frame()
  histvarquants_n <- 0
  for (histvar in unique(data$protein.name)) {
    histvardat <- data[data$protein.name == histvar, ]
    for (tissue in unique(data$cell.type...tissue)) {
      tisdat <- histvardat[histvardat$cell.type...tissue == tissue, ]
      for (timepoint in unique(data$timepoint)) {
        timedat <- tisdat[tisdat$timepoint == timepoint, ]
        for (repl in unique(data$biological.replicate)) {
          repldat <- timedat[timedat$biological.replicate == repl, ]
          histvarquant <- sum(repldat$quantification)
          histvarquants_n <- histvarquants_n + 1
          histvarquants[histvarquants_n,1] <- histvar
          histvarquants[histvarquants_n,2] <- tissue
          histvarquants[histvarquants_n,3] <- timepoint
          histvarquants[histvarquants_n,4] <- histvarquant
          histvarquants[histvarquants_n,5] <- repl
        }
      }
    }
  }
  colnames(histvarquants) <- c("protein.name", "cell.type...tissue", "timepoint", "histvarquant", "biological.replicate")
  return(histvarquants)
}

norm_by_histvar <- function(data) {
  # calculate histone variant abundances for all conditions
  histvarquants <- histvarquant(data)
  # normalize peptide abundances by histone variant abundances
  adj_data <- merge(data, histvarquants)
  adj_data$quantification <- adj_data$quantification / adj_data$histvarquant
  adj_data$histvarquant <- NULL
  return(adj_data)
}

sum_histvar_quants <- function(data) {
  # for each sample, sum up quantifications of both histone variants to get values for H3 total
  sumdata <- aggregate(data$quantification,
                       by = data[c("cell.type...tissue", "modifications", "timepoint", "biological.replicate")], sum)
  names(sumdata)[names(sumdata) == 'x'] <- "quantification"
  sumdata$protein.name <- "H3"
  return(sumdata)
}

prepPTMdata <- function(csv, histvars = TRUE, avrepls = TRUE) {
  # load data derived from CrosstalkDB (or in the respective format) and remove unnecessary columns
  # normalize according to histone variants, unless histvars = FALSE (for plotting H3 total data)
  # average, unless avrepls = FALSE (for comparison of replicates)
  
  data <- read.csv(csv, stringsAsFactors = FALSE)
  okdata <- cleanup_ctdb(data)
  okdata_0 <- zero_imp(okdata)
  
  if (histvars == TRUE) {
    adj_data <- norm_by_histvar(okdata_0)
  } else {
    adj_data <- sum_histvar_quants(okdata_0)
  }
  
  if (avrepls == TRUE) {
    adj_data <- av_repls(adj_data)
  } else {
    adj_data <- adj_data
  }
  
  return(adj_data)
  
}


###############################
## PTM abundance calculation ##
###############################

calcPTMab <- function(pepdata, outdir = getwd()) {
  
  # make list of all individual modifications present in the dataset
  modlist <- splitCombMod(pepdata$modifications)
  
  # calculate PTM individual and combinatorial abundances
  ab <- data.frame()
  ab_n <- 0
  for (hist in unique(pepdata$protein.name)) {
    histdat <- pepdata[pepdata$protein.name == hist, ]
    for (tissue in unique(pepdata$cell.type...tissue)) {
      tisdat <- histdat[histdat$cell.type...tissue == tissue, ]
      for (timepoint in unique(pepdata$timepoint)) {
        timedat <- tisdat[tisdat$timepoint == timepoint, ]
        for (repl in unique(pepdata$biological.replicate)) {
          repldat <- timedat[timedat$biological.replicate == repl, ]
          for (mi in modlist) {
            
            # get all combinatorial PTMs including that modification
            mi_combmod <- repldat[grep(mi, repldat$modifications),]
            
            # calculate relative abundance p_i of PTM m_i
            pi <- sum(mi_combmod$quantification)
            if (pi == 0) {
              next
            }

            # collect all interacting individual modifications m_j (this includes m_i)
            mjs_uniq <- splitCombMod(mi_combmod$modifications)

            for (mj in mjs_uniq) {
              
              if (mj == mi) {
                next
              }

              # calculate relative abundance p_j of PTM m_j
              pj <- sum(repldat[grep(mj, repldat$modifications),"quantification"])
              # calculate co-occurrence p_ij
              pij <- sum(mi_combmod[grep(mj, mi_combmod$modifications),"quantification"])
              if (pj == 0 | pij == 0) {
                next
              }
              # transform abundances
              pi_hat <- pi * (pi - pij) / ((1-pj) * pij)
              pj_hat <- pj * (pj - pij) / ((1-pi) * pij)
              # calculate interplay score from transformed abundances
              I <- I_trans(pi_hat, pj_hat)

              # add all values in a new row to ptm dataframe
              ab_n <- ab_n+1
              values <- list(hist, tissue, timepoint, repl, mi, pi, mj, pj, pij, pi_hat, pj_hat, I)
              for (colnr in seq(1, length(values))) {
                ab[ab_n,colnr] <- values[colnr]
              }
            }
          }
        }
      }
    }
  }
  
  colnames(ab) <- c("hist", "tissue", "timepoint", "repl", "mi", "pi", "mj", "pj", "pij", "pi_hat", "pj_hat", "I")
  
  # print data to file
  write.table(ab, paste0(outdir, "/ptm-abundances.tab"), quote = FALSE, sep = "\t", row.names = FALSE)
  
  return(ab)
  
}


##############
## Plotting ##
##############

raster <- function(start, end) {
  # generate data for rasterplot with transformed abundances and interplay score
  raster_df <- data.frame()
  n <- 0
  for (pi_t in log_seq(start, end, 100)) {
    for (pj_t in log_seq(start, end, 100)) {
      n <- n+1
      raster_df[n,1] <- pi_t
      raster_df[n,2] <- pj_t
      raster_df[n,3] <- I_trans(pi_t, pj_t)
    }
  }
  colnames(raster_df) <- c("pi_t", "pj_t", "I")
  return(raster_df)
}

encode <- function(data, splitplot_by = splitplot_by, connected = connected, group_by = group_by, colcode = colcode) {
  ## what should be encoded how? needs to be given as function argument in main function (CrossTalkMap())
  # splitplot_by: individual plots for which variable? (H3/histvars: tissue, repls: time; histone variants are taken care of automatically)
  # connected: the different instances of this variable are connected / grouped (H3/histvars: time, repls: replicates)
  # group_by: data points are grouped by this variable (H3/histvars: irrelevant (implicitely replicates), repls: tissue)
  # colcode: this variable is color-coded (H3/histvars: p_j, repls: tissue)
  data$splitplot_by <- data[[splitplot_by]]
  data$connected <- data[[connected]]
  data$group_by <- data[[group_by]]
  data$colcode <- data[[colcode]]
  return(data)
}

ptm_combs <- function(indmods, minlen = 2, maxlen = 2) {
  # function ptm_combs returns all possible combinations of a vector of individual ptms <indmods> up to a length of <maxlen> ptms
  # ptms within the combinatorial ptms are in order according to the peptide sequence
  # combinatorial ptms do not contain multiple ptms at the same position
  # default parameters produce all binary PTMs
  
  # transform to dataframe and add position
  modind_df <- as.data.frame(indmods)
  modind_df$pos <- as.integer(gsub("[[:upper:]]([[:digit:]]+)[[:alpha:]]+[[:digit:]]?$", "\\1", modind_df$indmods))
  modind_df <- modind_df[order(modind_df$pos),]
  
  # go through list, extract all other positions and concatenate each PTM to current PTM
  if (minlen == 1) {
    final_list <- indmods     # final list of combinations begins with all individual modifications
  } else {
    final_list <- vector()
  }
  iterate_over <- indmods
  samelen <- vector()
  nr_comb_ptms <- 2 # individual PTMs are already in the list (if minlen = 1), now start with adding binary PTMs
  
  while (nr_comb_ptms <= maxlen) {
    
    for (ptm in iterate_over) {
      
      # extract position of the PTM, match last sequence position
      position <- as.integer(regmatches(ptm, regexec("[[:upper:]]([[:digit:]]+)[[:alpha:]]+[[:digit:]]?$", ptm))[[1]][2])
      
      # get all PTMs at all downstream positions
      # (prevents generation of impossible combinations, i.e. different PTMs at the same position,
      # as well as redundant combinations, i.e. keeps the order)
      allothers <- modind_df[modind_df$pos > position, ]
      
      if (length(rownames(allothers)) == 0) {
        # last ptm has been in all combinations already, there is no modified position downstream
        next
      }
      
      # concatenate PTMs
      catptms <- paste0(ptm, allothers$indmods)
      samelen <- c(samelen, catptms)
      
    }
    
    # update list to iterate over
    iterate_over <- samelen
    # append to final list
    final_list <- c(final_list, samelen)
    # reset
    nr_comb_ptms <- nr_comb_ptms + 1
    samelen <- NULL
    
  }
  return(final_list)
}

make_point_labels <- function(data, which_label = which_label) {
  ## order data and generate labels only for the first data point in each subgroup
  data$timepoint <- as.numeric(gsub(" months", "", data$timepoint))
  data_ordered <- data[order(data$mi, data$mj, data$tissue, data$timepoint, data$repl),]
  data_labeled <- data.frame()
  for (mi in unique(data_ordered$mi)) {
    for (mj in unique(data_ordered$mj)) {
      for (group_by in unique(data_ordered$group_by)) {
        subgroup <- data_ordered[data_ordered$mi == mi & data_ordered$mj == mj & data_ordered$group_by == group_by,]
        if (nrow(subgroup) == 0) {
          next
        }
        subgroup$label <- ""
        if (which_label == "mj") {
          subgroup$label[1] <- subgroup$mj[1]
        } else if (which_label == "mimj") {
          mimj_ordered <- ptm_combs(c(mi, mj))
          subgroup$label[1] <- mimj_ordered
        }
        data_labeled <- rbind(data_labeled, subgroup)
      }
    }
  }
  return(data_labeled)
}

base_plot <- function(raster_df, start, end, hide_axes) {
  # generate base rasterplot with or without axes
  p <- ggplot(raster_df, aes(pi_t, pj_t)) +
    coord_cartesian(xlim = c(start, end), ylim = c(start, end)) +
    geom_raster(aes(fill = I)) +
    # line where interplay = 0
    geom_contour(aes(z = I), breaks = 0, size = 0.5, color = "gray70") +
    scale_x_continuous(trans = reverselog_trans(10), breaks = base_breaks(5), labels = prettyNum) +
    scale_y_continuous(trans = reverselog_trans(10), breaks = base_breaks(5), labels = prettyNum) +
    theme(axis.text=element_text(size=14, color = "black"),
          axis.title=element_text(size=14),
          plot.title = element_text(size = 18, hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5),
          legend.title = element_text(size = 14),
          legend.text = element_text(size = 14),
          panel.background = element_rect(fill = "white"),
          panel.border = element_rect(fill = NA, colour = "black"))
  if (hide_axes == TRUE) {
    p <- p + labs(x = "Transformed abundance PTM1", y = "Transformed abundance PTM2") +
      theme(axis.text = element_blank(), axis.ticks = element_blank())
  } else {
    p <- p + labs(x = expression(hat(p[i])), y = expression(hat(p[j])))
  }
  return(p)
}

add_title <- function(base_plot, mi, cond, hist, which_label) {
  # if data points are labeled with m_j only, m_i should be in the title
  # otherwise, labels are of the form mimj, then only condition and histone type in title
  if (which_label == "mj") {
    p <- base_plot + ggtitle(mi, subtitle = paste(cond, hist, sep = ", "))
  } else {
    p <- base_plot + ggtitle(paste(cond, hist, sep = ", "))
  }
  return(p)
}

add_interplay_col <- function(base_plot, col_scheme = "standard") {
  # determine color scheme for interplay score (plot background), where "standard" (default): gray gradient, "legacy": yellow to red
  if (col_scheme == "standard") {
    p <- base_plot + scale_fill_gradient2(low = "gray85", mid = "white", high = "gray85")
  } else if (col_scheme == "legacy") {
    p <- base_plot + scale_fill_gradient2(low = "#ffd27f", mid = "white", high = "#ff7f7f")
  } else {
    warning("Unknown argument to col_scheme")
    return()
  }
  return(p)
}

make_contour_label <- function(val, contour_label_form) {
  # define form of labels
  # default "short" only prints the value itself, "long" prepends with "I = "
  if (contour_label_form == "long") {
    contour_label <- paste("I =", val)
  } else {
    contour_label <- val
  }
  return(contour_label)
}

trans_ab2 <- function(trans_ab1, I) {
  # calculate transformed abundance 1 based on transformed abundance 2 and I
  pi <- (1 / (trans_ab1 * exp(1)^I))
  return(pi)
}

add_contour_labels <- function(plot, contour_breaks, raster_start, raster_end, contour_label_form) {
  # takes a list of contour line values and labels them at the bottom and right side of the plot
  for (I in contour_breaks) {
    # calculate x and y coordinates
    x_ab_trans <- max(trans_ab2(trans_ab1 = raster_end, I = I), raster_start)
    y_ab_trans <- min(trans_ab2(trans_ab1 = raster_start, I = I), raster_end)
    # make label (short or long form)
    contour_label <- make_contour_label(I, contour_label_form = contour_label_form)
    # add label to plot
    plot <- plot + annotate("text", label = contour_label, x = x_ab_trans, y = y_ab_trans, vjust = 0.5, hjust = 1, size = 3, angle = -45)
  }
  return(plot)
}

round2wards0 <- function(val, round2closest = 5) {
  # rounds value to the closest multiple of <round2closest> towards 0 (round up for negative and down for positive values)
  if (val < 0) {
    val_rounded <- ceiling(val / round2closest)*round2closest
  } else {
    val_rounded <- floor(val / round2closest)*round2closest
  }
  return(val_rounded)
}

add_contours <- function(base_plot, I_data, interval_size = 5, raster_start, raster_end, contour_label_form) {
  # add contour lines for interplay score at intervals of 5 (default)
  imin <- min(I_data$I)
  imax <- max(I_data$I)
  contour_breaks <- seq(round2wards0(imin, interval_size), round2wards0(imax, interval_size), interval_size)
  p <- base_plot + geom_contour(aes(z = I), show.legend = TRUE, size = 0.25, linetype = "longdash", color = "gray70",
                                breaks = contour_breaks)
  p <- add_contour_labels(p, contour_breaks, raster_start, raster_end, contour_label_form)
  return(p)
}

add_point_col <- function(base_plot, subgroup_data, all_data, colcode, col_scheme) {
  # determine which color scale to add for data points based on data type and background color scheme
  if (is.numeric(subgroup_data[,colcode]) == TRUE) {
    if (col_scheme == "standard") {
      p <- base_plot + scale_color_gradientn(colours = c("#00BFC4", "#3F5EFB", "#B06AB3", "#FC466B"),
                                             limits = c(min(all_data$colcode), max(all_data$colcode)))
    } else if (col_scheme == "legacy") {
      p <- base_plot + scale_color_gradientn(colours = c("#00BFC4", "#1A2980", "#B06AB3"),
                                             limits = c(min(all_data$colcode), max(all_data$colcode)))
    } else {
      warning("Unknown argument to col_scheme")
      return()
    }
    p <- p + guides(color = guide_colorbar(order = 1), fill = guide_colorbar(order = 2))
  } else {
    p <- base_plot + scale_color_discrete() +
      guides(color = guide_legend(order = 1), fill = guide_colorbar(order = 2))
  }
  return(p)
}

add_col_legend_label <- function(base_plot, colcode) {
  # formatting of color-code legend label
  # assuming colcode is "pxyz", returns labels as "p[xyz]" in the plot legend (xyz as subscript to p)
  subscr <- gsub("p", "", colcode)
  col_label <- bquote(p[.(subscr)]) #expression(p[subscr])
  p <- base_plot + labs(color = col_label)
  return(p)
}

add_points <- function(base_plot, data) {
  # add data points and labels to plot
  p <- base_plot +
    geom_point(data = data, mapping = aes(x = pi_hat, y = pj_hat, color = colcode), size = 1) +
    geom_text_repel(data = data, mapping = aes(x = pi_hat, y = pj_hat, label = label, color = colcode),
                    size = 3, segment.size = 0.025)
  return(p)
}

add_paths <- function(plot, data, with_arrows) {
  # add paths beween data points, if specified with arrows
  for (mi in sort(unique(data$mi))) {
    for (mj in sort(unique(data$mj))) {
      mimj_connected <- data[data$mi == mi & data$mj == mj,]
      # if mj present at at least two time points, connect them by arrow (or for replicates by line only)
      if (nrow(mimj_connected) > 1 && with_arrows == TRUE) {
        plot <- plot + geom_path(data = mimj_connected, mapping = aes(x = pi_hat, y = pj_hat, color = colcode), size = 0.5,
                                   arrow = arrow(length = unit(0.18, "cm")))
      } else if (nrow(mimj_connected) > 1 && with_arrows == FALSE) {
        plot <- plot + geom_path(data = mimj_connected, mapping = aes(x = pi_hat, y = pj_hat, color = colcode), size = 0.5)
      }
    }
  }
  return(plot)
}

plot_all <- function(plotlist_all, ptm_data, outdir, splitplot_by, filename_string, filename_ext) {
  ## generate multiplot from list of plots, including .tab file of contained data
  ## ptm_data is used to calculate the theoretical number of individual plots, since some in the list might be empty
  
  # calculate number of multi plot columns and rows
  # there will be one plot for each tissue and each histone (variant)
  nr_cond <- length(unique(ptm_data$splitplot_by))
  nr_histvar <- length(unique(ptm_data$hist))
  plotsn <- nr_cond * nr_histvar
  nrrows <- floor(sqrt(plotsn))
  nrcols <- ceiling(plotsn / nrrows)
  
  # define layout for multiplot  
  if (nr_histvar == 1) {
    # only one histone variant, disregard this dimension
    lay <- vector()
    for (rownr in seq(1, nrrows)) {
      plot_last <- nrcols * rownr
      plot_1st <- plot_last - nrcols + 1
      row <- seq(plot_1st, plot_last)
      lay <- rbind(lay, row)
    }
  } else {
    # one row per histone variant, one column per tissue/condition
    lay <- vector()
    for (colnr in seq(1, nrcols)) {
      plot_last <- nrrows * colnr
      plot_1st <- plot_last - nrrows + 1
      col <- seq(plot_1st, plot_last)
      lay <- cbind(lay, col)
    }
  }
  p <- grid.arrange(grobs = plotlist_all, ncol = nrcols, nrow = nrrows, layout_matrix = lay)
  
  if (is.null(filename_string)) {
    # return plot object
    return(p)
  } else {
    # determine histone type for file name
    if (nrow(ptm_data[grep('\\.', ptm_data$hist),]) > 0) {
      filename_hist <- "histvars"
    } else {
      filename_hist <- "H3"
    }
    
    # calculate width and height of multi plot
    plotunit <- 6 # each individual plot should be 6x6 cm
    plotheight <- nrrows * plotunit
    plotwidth <- nrcols * plotunit
    
    # plot
    ggsave(filename = paste0(outdir, "/crosstalkmap_splitby-", splitplot_by, "_", filename_hist, "_", filename_string, ".", filename_ext), plot = p,
        width = plotwidth, height = plotheight)

    # print values for all plotted PTMs to file
    ptm_data_sort <- ptm_data[with(ptm_data, order(hist, tissue, mj, timepoint)),]
    write.table(ptm_data_sort,
                paste0(outdir, "/crosstalkmap_splitby-", splitplot_by, "_", filename_hist, "_", filename_string, ".tab"),
                quote = FALSE, sep = "\t", row.names = FALSE)
  }
}

CrossTalkMap <- function(ptm_data, splitplot_by = "tissue", colcode = "pj", connected = "timepoint", group_by = "repl",
                         connect_dots = TRUE, with_arrows = TRUE, which_label = "mj",
                         col_scheme = "standard", contour_lines = TRUE, contour_labels = c("short", "long"), hide_axes = TRUE,
                         filename_string = NULL, filename_ext = "pdf", outdir = getwd()) {
  ## take data frame with transformed PTM abundances as input
  ## crosstalk maps for all PTM combinations of m_i and m_j contained in the input data frame are plotted
  ## encoding according to splitplot_by, connected, group_by, colcode (see encode())
  ## which_label defines labels for groups of data points: "mj" (default) for individual or "mimj" for combinatorial PTM
  ## col_scheme determines the color of the interplay score gradient (plot background), where "standard" (default): gray scale, "legacy": yellow to red
  ## by default, returns plot object
  ## otherwise, if argument to filename_string is provided, a file name is constructed and a file created
  ## (default: pdf, otherwise one of "eps", "ps", "tex" (pictex), "jpeg", "tiff", "png", "bmp", "svg" or "wmf", might require additional packages)
  ## filename_string will be attached to the plot file name starting with crosstalkmap_<mi> and to reflect how m_js to plot have been selected
  
  # define start and end for x and y scales based on min and max in pi_hat and pj_hat
  # will be the same for all plots to make them comparable
  # start cannot be 0 (this would be the case if pi = pij); therefore, find smallest value > 0 for min
  start <- min(ptm_data$pi_hat[is.finite(ptm_data$pi_hat) & ptm_data$pi_hat > 0],
               ptm_data$pj_hat[is.finite(ptm_data$pj_hat) & ptm_data$pj_hat > 0])
  end <- max(ptm_data$pi_hat[is.finite(ptm_data$pi_hat)], ptm_data$pj_hat[is.finite(ptm_data$pj_hat)])
  
  # generate data for rasterplot
  raster_df <- raster(start, end)
  
  # what should be encoded how? copy respective columns
  ptm_data <- encode(ptm_data, splitplot_by = splitplot_by, connected = connected, group_by = group_by, colcode = colcode)

  # print all x plots (histone (variants) * 4 tissues) on one page
  plotlist_all <- list()
  plot_count_all <- 0

  # iterate over individual conditions, making one plot for each
  for (cond in unique(ptm_data[[splitplot_by]])) {
    split_plot <- ptm_data[ptm_data[[splitplot_by]] == cond, ]
    for (hist in unique(ptm_data$hist)) {

      hist_plot <- split_plot[split_plot$hist == hist, ]
      plot_count_all <- plot_count_all + 1

      # make labels and order data accordingly
      hist_labeled <- make_point_labels(hist_plot, which_label = which_label)
      # make raster plot
      p <- base_plot(raster_df, start, end, hide_axes = hide_axes)
      # add title
      p <- add_title(p, mi, cond, hist, which_label)
      # add color scale for interplay score / raster plot
      p <- add_interplay_col(p, col_scheme)
      # add interplay score contour lines
      if (contour_lines == TRUE) {
        contour_label_form <- match.arg(contour_labels)
        p <- add_contours(p, raster_df, raster_start = start, raster_end = end, contour_label_form = contour_label_form)
      }
      # determine which color scale to add for data points
      p <- add_point_col(p, subgroup_data = hist_labeled, ptm_data, colcode, col_scheme)
      # format color-code legend label
      p <- add_col_legend_label(p, colcode)
      # add data points and labels
      p <- add_points(p, hist_labeled)
      # add paths between points
      if (connect_dots == TRUE) {
        p <- add_paths(p, hist_labeled, with_arrows)
      }

      # add plot to list
      plotlist_all[[plot_count_all]] <- p

    }
  }

  # generate multiplot from list of plots, incl .tab of all data points contained in the plot
  plot_all(plotlist_all, ptm_data, outdir, splitplot_by, filename_string, filename_ext)

}

################
## Line plots ##
################

## For the standard case (over time in each tissue)

line_ab <- function(data, outdir = getwd()) {
  ## line plot for abundance change over time
  ## layout for manual integration into crosstalk map, hence no tissue / histone variant labels
  data$timepoint_nr <- as.numeric(gsub(" months", "", data$timepoint))
  data <- data[order(data$timepoint_nr),]
  tissue <- unique(data$tissue)
  hist <- unique(data$hist)
  if (length(tissue) > 1 || length(hist) > 1) {
    warning("Provide data for only one tissue and histone (variant)!")
    return()
  }
  # plot
  pdf(paste0(outdir, "/ab_", unique(data$mi), "_", tissue, "_", sub(" ", "-", hist), ".pdf"), height = 4.7)
  p <- ggplot(data, aes(timepoint_nr, pi)) +
    geom_line(size = 1.5) +
    geom_point(size = 4) +
    labs(x = "Months", y = paste(unique(data$mi), "abundance")) +
    scale_x_continuous(breaks = data$timepoint_nr, labels = data$timepoint_nr) +
    # to prevent some y tick labels to be cut off (although it still cuts off in one or the other case)
    scale_y_continuous(limits = c(min(data$pi), max(data$pi) + ((max(data$pi)-min(data$pi))/100))) +
    theme_classic() +
    theme(text = element_text(size = 30)) +
    coord_cartesian(clip = "off")
  print(p)
  dev.off()
}

line_ct <- function(data, outdir = getwd()) {
  ## line plot for changes of abundances, co-occurrences and interplay score over time for PTM pair in one tissue
  ## if no defined Interplay score at an individual time point, all other measures are plotted,
  ## but the interplay score data point is missing
  data$timepoint_nr <- as.numeric(gsub(" months", "", data$timepoint))
  data <- data[order(data$timepoint_nr),]
  tissue <- unique(data$tissue)
  hist <- unique(data$hist)
  if (length(tissue) > 1 || length(hist) > 1) {
    warning("Provide data for only one tissue and histone (variant)!")
    return()
  }
  # plot
  pdf(paste0(outdir, "/ct-params_", data$mi, data$mj, "_", tissue, "_", sub(" ", "-", hist), ".pdf"))
  par(xpd = TRUE, mar = c(6,5,2,5))
  plot(data$timepoint_nr, data$pi, type = "b", col = "orange", lwd = 2,
       xlim = c(min(data$timepoint_nr), max(data$timepoint_nr)), ylim = c(0, max(c(data$pi, data$pj))),
       xlab = "Months", ylab = "Abundances", xaxt = "n")
  title(paste0(toupper(substr(tissue, 1, 1)), substr(tissue, 2, nchar(tissue)), ", ", hist), line = 0.4)
  axis(side = 1, at = data$timepoint_nr)
  lines(data$timepoint_nr, data$pj, type = "b", col = "blue", lwd = 2)
  lines(data$timepoint_nr, data$pij, type = "b", col = "red", lwd = 2)
  # if interplay score is not defined in any sample (i.e. pi_hat or pj_hat = 0 in all samples),
  # then pi or pj == pij in all samples, and no defined axis limits for score
  # in that case, don't plot score axis and indicate which abundance equals co-occurrence
  if (all(data$pi_hat == 0) || all(data$pj_hat == 0)) {
    if (all(data$pi_hat == 0)) {
      equal_pij <- unique(data$mi)
    } else if (all(data$pj_hat == 0)) {
      equal_pij <- unique(data$mj)
    }
    mtext(side = 3, line = -1.2, adj = 0.95, paste0("freq(", equal_pij, ") == freq(", unique(data$mi), unique(data$mj), ")"))
  } else {
    par(new = T)
    plot(data$timepoint_nr, data$I, type = "b", col = "black", lwd = 2,
         xlim = c(min(data$timepoint_nr), max(data$timepoint_nr)),
         axes = F, xlab = NA, ylab = NA)
    axis(side = 4)
    mtext(side = 4, line = 3, "Interplay Score")
  }
  legend("bottom", inset = c(0, -0.225), horiz = TRUE, lty = "solid",
         col = c("orange", "blue", "red", "black"), lwd = 2,
         legend = c(data$mi[1], data$mj[1], paste0(data$mi[1], data$mj[1]), "Interplay"), bty = "n")
  dev.off()
}


###############
## Utilities ##
###############

splitCombMod <- function(comblist){
  # comblist is a vector containing combinatorial modifications
  # combinatorial modifications will be converted into a unique list of all contained individual modifications
  mod_ind <- gsub("(?!^)(?=[[:upper:]])", " ", comblist, perl=T) # split into individual modifications (https://stackoverflow.com/questions/7988959/splitting-string-based-on-letters-case)
  mod_uniq <- unique(unlist(strsplit(mod_ind, " "))) # separate and make unique list
  return(mod_uniq)
}

I_trans <- function(p_i_t, p_j_t, epsilon = 0.00000000001){
  # calculates the interplay score from transformed individual PTM abundances p_i_t and p_j_t
  if(abs(p_i_t) < epsilon || abs(p_j_t) < epsilon){
    return(NaN)
  } else {
    return(log(1 / (p_i_t * p_j_t)))
  }
}

log_seq <- function(lbound, ubound, stepnr) {
  # generates a sequence equally spaced on logarithmic scale from lower bound to upper bound with x steps
  if ((ubound / lbound) > 0 && stepnr > 0) {
    k <- (ubound / lbound)^(1/stepnr)
    l <- log(lbound, k)
    logseq <- sapply(seq(0, stepnr), function(x) (k^(l + x)))
    return(logseq)
  }
}

base_breaks <- function(n = 10){
  # function for nice axis breaks (https://stackoverflow.com/a/22227846)
  # argument determines approximate number of labeled tick marks, default = 10
  function(x) {
    axisTicks(log10(range(x, na.rm = TRUE)), log = TRUE, n = n)
  }
}
