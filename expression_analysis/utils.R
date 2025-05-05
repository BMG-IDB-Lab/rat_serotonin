`%not in%` <- Negate(`%in%`)

myFeaturePlot <- function(object, features, 
                          assay = DefaultAssay(object), slot = "data", 
                          ncol = NULL, 
                          title.size = 10,
                          keep.scale.all = F,
                          title.face = "bold",
                          keep.axis.titles = F,
                          keep.scale.bar = T,
                          keep.axis.numbers = T,
                          keep.axis.ticks = T,
                          keep.axis.lines = T,
                          custom.cols = NULL,
                          legend.key.size = NULL,
                          ...) {
  
  if (is.null(ncol)) {ncol <- max(1, length(features) %/% 2)}
  DefaultAssay(object) <- assay
  p <- FeaturePlot(object, slot = slot, 
                   features = features, 
                   combine = FALSE,
                   ...)
  
  for(i in 1:length(p)) {
    p[[i]] <- p[[i]] + theme(title = element_text(size = title.size),
                             panel.border = element_blank())
  }
  
  if (!keep.axis.titles) {
    for(i in 1:length(p)) {
      p[[i]] <- p[[i]] + theme(axis.title.x = element_blank(),
                               axis.title.y = element_blank(),
                               legend.key.size = unit(0.5, 'cm'))
    }
  }
  if (keep.scale.all) {
    for(i in 1:length(p)) {
      p[[i]] <- p[[i]] + # theme(legend.key.size = unit(0.5, 'cm')) +
        scale_color_gradient(low = "lightgrey", high = "blue", 
                             limits = c(0,max(GetAssayData(object, 
                                                           assay = assay, 
                                                           slot = slot)[features, ])))
    }
  }
  if (!keep.scale.bar) {
    for (i in 1:length(p)) {
      p[[i]] <- p[[i]] + guides(color = "none")
    }
  }
  if (!keep.axis.numbers) {
    for (i in 1:length(p)) {
      p[[i]] <- p[[i]] + theme(axis.text.x = element_blank(),
                               axis.text.y = element_blank())
    }
  }
  if (!keep.axis.ticks) {
    for (i in 1:length(p)) {
      p[[i]] <- p[[i]] + theme(axis.ticks = element_blank())
    }
  }
  if (!keep.axis.lines) {
    for (i in 1:length(p)) {
      p[[i]] <- p[[i]] + theme(axis.line = element_blank())
    }
  }
  for (i in 1:length(p)) {
    p[[i]] <- p[[i]] + theme(plot.title = element_text(face = title.face))
  }
  if (length(custom.cols) > 0) {
    for (i in 1:length(p)) {
      p[[i]] <- p[[i]] + scale_colour_gradientn(colors = custom.cols)
    }
  }
  if(!is.null(legend.key.size)) {
    for (i in 1:length(p)) {
      p[[i]] <- p[[i]] + theme(legend.key.size = legend.key.size)
    }
  }
  
  cowplot::plot_grid(plotlist = p, ncol = ncol)
  
}


proportionsPlot <- function(obj, 
                            condition = "treatment",
                            samples = "sample",
                            cell_type = "Cell_type",
                            condition1, condition2,
                            sample_levels, ct_levels,
                            group_color_4var = F, 
                            ngroups1 = c(2,2), ngroups2 = c(2,2),
                            rev_group_color = F) {
  meta <- obj@meta.data[, c(condition, samples, cell_type)]
  meta$cluster <- meta[, cell_type]
  meta$treatment <- meta[, condition]
  meta$sample <- meta[, samples]
  meta$cell <- rownames(meta)
  
  meta <- dplyr::select(meta, cluster, treatment, sample, cell)
  
  meta.pct <- meta %>% group_by(treatment, sample) %>%
    mutate(n = n()) %>%
    ungroup() %>%
    group_by(cluster, sample) %>%
    summarise(pct = (n() / n) * 100,
              treatment) %>%
    distinct()
  
  
  p1 <- ggplot(meta.pct, aes(x = cluster, y = pct, fill = sample, width = 0.7)) +
    geom_bar(stat = "identity", position = position_dodge(preserve = "single")) +
    facet_wrap(. ~ treatment) +
    theme_bw() + RotatedAxis() +
    ggtitle(sprintf("Percentage of cells in %s/%s", condition1, condition2))
  
  
  meta.pct$sample <- factor(meta.pct$sample, levels = sample_levels)
  meta.pct$cluster <- factor(meta.pct$cluster, levels = ct_levels)
  
  if (group_color_4var) {
    meta.pct.groups <- meta %>% dplyr::select(sample, treatment) %>% distinct()
    rownames(meta.pct.groups) <- NULL
    # meta.pct <- meta.pct %>% tidyr::complete(sample,cluster)
    meta.pct <- merge(meta.pct, meta.pct.groups, by = "sample", all.y = T)
    meta.pct <- meta.pct %>% dplyr::select(-treatment.x) %>% dplyr::rename(treatment = treatment.y)
    meta.pct[is.na(meta.pct)] <- 0
    if (rev_group_color) {
      colors <- c(rep(c("steelblue1", "dodgerblue3"), times = ngroups1),
                  rep(c("pink", "tomato"), times = ngroups2))
    } else {
      colors <- c(rep(c("pink", "tomato"), times = ngroups1), 
                  rep(c("steelblue1", "dodgerblue3"), times = ngroups2))
    }
    p2 <- ggplot(meta.pct, aes(x = cluster, y = pct, group = sample, fill = sample, width = 0.7)) +
      geom_bar(stat = "identity", position = position_dodge(preserve = "single"), color = "white") +
      theme_bw() + RotatedAxis() +
      ggtitle(sprintf("Percentage of cells in %s/%s", condition1, condition2)) +
      ylab("Percent of cells from sample") +
      theme(axis.text.x = element_text(size = 12),
            axis.title.x = element_blank(),
            legend.position="bottom") + 
      guides(fill=guide_legend(title="Group")) +
      scale_fill_manual(values = colors)
  } else {
    print("p1")
    
    if (rev_group_color) {
      colors <- c("dodgerblue3", "tomato")
    } else {
      colors <- c("tomato", "dodgerblue3")
    }
    print(colors)
    
    meta.pct$sample <- factor(meta.pct$sample, levels = sample_levels)
    meta.pct$cluster <- factor(meta.pct$cluster, levels = ct_levels)
    
    p2 <- ggplot(meta.pct, aes(x = cluster, y = pct, group = sample, fill = treatment, width = 0.7)) +
      geom_bar(stat = "identity", position = position_dodge(preserve = "single"), color = "white") +
      theme_bw() + RotatedAxis() +
      ggtitle(sprintf("Percentage of cells in %s/%s", condition1, condition2)) +
      ylab("Percent of cells from sample") +
      theme(axis.text.x = element_text(size = 12),
            axis.title.x = element_blank(),
            legend.position="bottom") + 
      guides(fill=guide_legend(title="Group")) +
      scale_fill_manual(values = colors)
  }
  
  return(list(
    cond_split = p1,
    cond_sample = p2
  ))
  
}


cao_compositional_rmCluster <- function(cao, data, cluster) {
  
  if (TRUE) {
    type=NULL; name='cell.density'; size=1.5; palette=NULL;
    adjust.pvalues=NULL; contours=NULL; contour.color='black'; contour.conf='10%';
    plot.na=FALSE; color.range=NULL; mid.color='gray95';
    scale.z.palette=adjust.pvalues; min.z=qnorm(0.9)
    
    palette <- c('blue', mid.color, 'red')
    palette %<>% grDevices::colorRampPalette(space="Lab")
    
    dens.res <- cao[["test.results"]][["cell.density"]]
    type <- if (length(dens.res$diff) == 0) 'permutation' else names(dens.res$diff)[1]
    
    scores <- dens.res$diff[[type]]
    if (is.null(adjust.pvalues)) {
      if (!is.null(scores$adj)) {
        scores <- scores$adj
        adjust.pvalues <- TRUE
      } else {
        scores <- scores$raw
        adjust.pvalues <- FALSE
      }
    } else if (adjust.pvalues) {
      if (is.null(scores$adj)) {
        warning("Adjusted scores are not estimated. Using raw scores. ",
                "Please, run estimateCellDensity with adjust.pvalues=TRUE")
        scores <- scores$raw
      } else {
        scores <- scores$adj
      }
    } else {
      scores <- scores$raw
    }
    
    if (dens.res$method == 'graph') {
      density.emb <- cao$embedding
      scores %<>% .[intersect(names(.), rownames(density.emb))]
    } else if (dens.res$method == 'kde') {
      density.emb <- dens.res$density.emb[,1:2]
    } else stop("Unknown method: ", dens.res$method)
    
    # remove cluster
    cells_rm <- Idents(data)
    cells_rm <- names(cells_rm[which(cells_rm == cluster)])
    scores <- scores[which(names(scores) %not in% cells_rm)]
    orig.emb <- cao[["embedding"]]
    cao[["embedding"]] <- cao[["embedding"]][setdiff(rownames(cao[["embedding"]]),
                                                     cells_rm),]

    scores <- sort(scores)
    density.mat <- dens.res$density.mat[names(scores),]
    density.emb <- density.emb[names(scores),]
    
    if (is.null(color.range)) {
      color.range <- c(-1, 1) * max(abs(scores))
    } else {
      color.range %<>% cacoa:::parseLimitRange(scores)
      color.range <- max(abs(color.range)) %>% {c(-., .)}
      scores %<>% pmin(color.range[2]) %>% pmax(color.range[1])
    }
    
    leg.title <- if (type == 'subtract') 'Prop. change' else {if (adjust.pvalues) 'Z adj.' else 'Z-score'}
    gg <- cao$plotEmbedding(density.emb, colors=scores, size=size, legend.title=leg.title, palette=palette,
                            midpoint=0, plot.na=plot.na, color.range=color.range)
    
    gg$scales$scales %<>% .[sapply(., function(s) !("colour" %in% s$aesthetics))]
    gg <- gg + cacoa:::getScaledZGradient(min.z=min.z, palette=palette, color.range=color.range)
    p1 <- gg + theme_void()
    
    gg$layers[[1]]$aes_params$size <- 3.5
    p2 <- gg + theme_void() +
      theme(legend.key.size = unit(1, "cm"),
            legend.text = element_text(size = 20),
            legend.title = element_text(size = 25))
    
    cao[["embedding"]] <- orig.emb
    
    return(list(plot = p1,
                plot_zoom = p2))
  }
  
}


cao_compositional_orderSize <- function(cao) {
  
  type=NULL; name='cell.density'; size=1.5; palette=NULL;
  adjust.pvalues=NULL; contours=NULL; contour.color='black'; contour.conf='10%';
  plot.na=FALSE; color.range=NULL; mid.color='gray95';
  scale.z.palette=adjust.pvalues; min.z=qnorm(0.9)
  
  palette <- c('blue', mid.color, 'red')
  palette %<>% grDevices::colorRampPalette(space="Lab")
  
  dens.res <- cao[["test.results"]][["cell.density"]]
  type <- if (length(dens.res$diff) == 0) 'permutation' else names(dens.res$diff)[1]
  
  scores <- dens.res$diff[[type]]
  if (is.null(adjust.pvalues)) {
    if (!is.null(scores$adj)) {
      scores <- scores$adj
      adjust.pvalues <- TRUE
    } else {
      scores <- scores$raw
      adjust.pvalues <- FALSE
    }
  } else if (adjust.pvalues) {
    if (is.null(scores$adj)) {
      warning("Adjusted scores are not estimated. Using raw scores. ",
              "Please, run estimateCellDensity with adjust.pvalues=TRUE")
      scores <- scores$raw
    } else {
      scores <- scores$adj
    }
  } else {
    scores <- scores$raw
  }
  
  if (dens.res$method == 'graph') {
    density.emb <- cao$embedding
    scores %<>% .[intersect(names(.), rownames(density.emb))]
  } else if (dens.res$method == 'kde') {
    density.emb <- dens.res$density.emb[,1:2]
  } else stop("Unknown method: ", dens.res$method)
  
  # sort
  scores <- c(scores[scores == 0],
              sort(scores[scores<0], decreasing = T),
              sort(scores[scores>0]))
  density.mat <- dens.res$density.mat[names(scores),]
  density.emb <- density.emb[names(scores),]
  
  if (is.null(color.range)) {
    color.range <- c(-1, 1) * max(abs(scores))
  } else {
    color.range %<>% cacoa:::parseLimitRange(scores)
    color.range <- max(abs(color.range)) %>% {c(-., .)}
    scores %<>% pmin(color.range[2]) %>% pmax(color.range[1])
  }
  
  leg.title <- if (type == 'subtract') 'Prop. change' else {if (adjust.pvalues) 'Z adj.' else 'Z-score'}
  gg <- cao$plotEmbedding(density.emb, colors=scores, size=size, legend.title=leg.title, palette=palette,
                          midpoint=0, plot.na=plot.na, color.range=color.range)
  
  gg$scales$scales %<>% .[sapply(., function(s) !("colour" %in% s$aesthetics))]
  gg <- gg + cacoa:::getScaledZGradient(min.z=min.z, palette=palette, color.range=color.range)
  p1 <- gg + theme_void()
  
  gg$layers[[1]]$aes_params$size <- 4
  p2 <- gg + theme_void() +
    theme(legend.key.size = unit(1, "cm"),
          legend.text = element_text(size = 20),
          legend.title = element_text(size = 25))
  
  return(list(p1 = p1, p2 = p2))
  
}


cao_expression_top <- function(cao, data,
                               top = 10,
                               cluster_rm = "51_Dop_Chol_Eef1b2.Cd59b") {
  
  shifts <- cao[["test.results"]][["cluster.free.expr.shifts"]]
  z.scores <- shifts$z_scores %>% sort(decreasing = F)
  shifts <- shifts$shifts_smoothed %>% sort(decreasing = F)
  
  color.range=c("0", "97.5%")
  shifts %<>% na.omit()
  color.range %<>% cacoa:::parseLimitRange(shifts)
  shifts %<>% pmax(color.range[1]) %>% pmin(color.range[2])
  
  # barplots
  df <- Reduce(function(x,y) merge(x,y,by="cell"), 
               list(as.data.frame(shifts) %>% tibble::rownames_to_column("cell"),
                    as.data.frame(z.scores) %>% tibble::rownames_to_column("cell"),
                    as.data.frame(data$res.combined_det_markers) %>% 
                      tibble::rownames_to_column("cell") %>% 
                      `colnames<-`(c("cell", "res")))) %>% 
    select(-cell) %>% 
    group_by(res) %>% 
    summarise(shift_mean = mean(shifts),
              z_mean = mean(z.scores))
  
  
  df_z <- df %>% slice_max(order_by = z_mean, n = top) %>% 
    filter(res != cluster_rm) %>% 
    arrange(z_mean) %>% as.data.frame() %>% 
    mutate(res = as.character(res),
           res = factor(res, levels = unique(res)))
  palette=cacoa::brewerPalette("YlOrRd", rev=FALSE)
  max.score=max(df_z$z_mean, na.rm=TRUE); min.z=min(df_z$z_mean)
  col.vals <- c(0, seq(min.z, max.score, length.out=20))
  max.score <- c(0, max.score)
  col.vals %<>% scales::rescale()
  
  p_barplot <- ggplot(df_z,
         aes(x = res, y = z_mean, fill = z_mean)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    theme_minimal() +
    scale_fill_gradientn(colors=palette(21), 
                         values=col.vals, 
                         limits=max.score, breaks = c(0,1,2)) +
    scale_y_continuous(breaks = c(0, 1, 2)) +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          
          axis.text.y = element_text(size = 20, color = "black"),
          axis.ticks.y = element_blank(),
          
          axis.text.x = element_text(size = 18, color = "black"),
          axis.line.x = element_line(),
          axis.ticks.x = element_line(),
          
          panel.grid = element_blank(),
          
          legend.position = "right",
          legend.text = element_text(size = 15),
          legend.title = element_text(size = 18, vjust = 0.8)) +
    theme(plot.margin = margin(0,70,0,0))
  
  
  # expression umap z-score
  emb_shift <- cao$embedding[names(shifts),]
  emb_z <- cao$embedding[names(z.scores),]
  
  palette=cacoa::brewerPalette("YlOrRd", rev=FALSE)
  ggs <- mapply(function(cls, lt, emb) {
    cao$plotEmbedding(emb, size = 1.2,
                      colors=cls, plot.na=FALSE, alpha=0.2, 
                      palette=palette, 
                      legend.title=lt) +
      theme(legend.background = element_blank())
  }, list(shifts, z.scores), c("Distance", "Z-score"), list(emb_shift, emb_z), 
  SIMPLIFY=FALSE)
  
  scale.z.palette=TRUE
  if (scale.z.palette) {
    ggs[[2]]$scales$scales %<>% .[sapply(., function(s) !("colour" %in% s$aesthetics))]
    max.score <- max(z.scores, na.rm=TRUE)
    ggs[[2]] <- ggs[[2]] + cacoa:::getScaledZGradient(min.z=qnorm(0.9), 
                                                      palette=palette, 
                                                      color.range=max.score)
  }
  
  ggs <- lapply(ggs, function(x) {
    x <- x + theme_void() +
      theme(plot.margin = margin(10,15,10,10),
            legend.key.size = unit(0.5, "cm"))
  })
  
  cell.groups=cao$cell.groups
  n <- names(cell.groups)
  cell.groups <- case_when(
    cell.groups %in% gsub("(\\d+)_[GABA|Glut].*", "\\1", df_z$res) ~ cell.groups,
    TRUE ~ ""
  )
  cell.groups <- setNames(cell.groups, n)
  
  font.size=c(9, 6)
  if (!is.null(cell.groups)) {
    ggs %<>% lapply(cacoa:::transferLabelLayer, 
                    cao$plotEmbedding(groups=cell.groups, size=1.2), 
                    font.size=font.size)
  }
  
  p_expr <- cowplot::plot_grid(plotlist=mapply(function(x, lab) {
    x$layers[[1]]$aes_params$size <- 2
    x <- x + theme(legend.key.size = unit(0.8, "cm"),
                   legend.text = element_text(size = 18),
                   legend.title = element_text(size = 22)) +
      theme(plot.margin = margin(10,5,5,5),
            legend.position = "bottom")
    cowplot::plot_grid(NULL, x, ncol = 1, rel_heights = c(0.08, 1),
                       labels = c(lab, NULL), label_size = 25)
  }, ggs, c("Shifts", "Adj. z-scores"), SIMPLIFY = F), ncol=2)
  
  return(list(plot_barplot = p_barplot,
              plot_zscore = p_expr))
  
}

cao_programs <- function(cao) {
  
  gene.progs <- cao[["test.results"]][["gene.programs"]]
  prog.ids <- (sapply(gene.progs$genes.per.clust, length) > 10) %>% which()
  
  plots_list <- list()
  for (i in seq_along(prog.ids)) {
    cur.scores <- gene.progs$program.scores[i,]
    cells_missing <- setdiff(Cells(data), names(cur.scores))
    scores_missing <- setNames(rep(0, length(cells_missing)),
                               cells_missing)
    cur.scores <- c(cur.scores, scores_missing)
    cur.range <- cacoa:::parseLimitRange(c("0.5%", "99.5%"), cur.scores)
    cur.scores %<>% pmax(cur.range[1]) %>% pmin(cur.range[2])
    
    pr_name <- paste0("Program ", i)
    df_plot <- merge(cao$embedding %>% as.data.frame() %>% 
                       tibble::rownames_to_column("cell"),
                     data.frame(cell = names(cur.scores),
                                score = cur.scores),
                     by = "cell")
    
    if (min(df_plot$score) < 0) {
      df_plot <- df_plot %>% arrange(-score)
      breaks <- quantile(scales::extended_breaks()(df_plot$score), probs = c(0.1, 0.5, 1))
      p <- ggplot(df_plot, aes(x = UMAP_1, y = UMAP_2, color = score)) +
        geom_point(size = 3, shape = 20, alpha = 0.5) +
        scale_color_gradient(low = "blue", high = "grey90", 
                             breaks =  breaks,
                             labels = breaks) +
        theme_void() +
        theme(plot.margin = margin(5,10,5,10)) +
        theme(legend.title=element_text(size=18, vjust = 1), 
              legend.text=element_text(size=15)) +
        theme(plot.title = element_text(size = 20,  
                                        hjust = 0.5)) +
        theme(legend.key.size = unit(0.3, "cm")) +
        ggtitle(pr_name)
    } else {
      df_plot <- df_plot %>% arrange(score)
      breaks <- quantile(scales::extended_breaks()(df_plot$score), probs = c(0, 0.5, 0.9))
      p <- ggplot(df_plot, aes(x = UMAP_1, y = UMAP_2, color = score)) +
        geom_point(size = 3, shape = 20, alpha = 0.5) +
        scale_color_gradient(low = "grey90", high = "red", 
                             breaks =  breaks,
                             labels = breaks) +
        theme_void() +
        theme(plot.margin = margin(5,10,5,10)) +
        theme(legend.title=element_text(size=18, vjust = 1), 
              legend.text=element_text(size=15)) +
        theme(plot.title = element_text(size = 20, 
                                        hjust = 0.5)) +
        theme(legend.key.size = unit(0.3, "cm")) +
        ggtitle(pr_name)
    }
    
    plots_list[[i]] <- p
    
  }
  
  return(plots_list)
  
}

cao_programs_barplot <- function(cao, data, anno = "res.combined_det_markers") {
  
  gene.progs <- cao[["test.results"]][["gene.programs"]]
  prog.ids <- (sapply(gene.progs$genes.per.clust, length) > 10) %>% which()
  
  plots_list <- list()
  for (i in seq_along(prog.ids)) {
    cur.scores <- gene.progs$program.scores[i,]
    cells_missing <- setdiff(Cells(data), names(cur.scores))
    scores_missing <- setNames(rep(0, length(cells_missing)),
                               cells_missing)
    cur.scores <- c(cur.scores, scores_missing)
    cur.range <- cacoa:::parseLimitRange(c("0.5%", "99.5%"), cur.scores)
    cur.scores %<>% pmax(cur.range[1]) %>% pmin(cur.range[2])
    
    pr_name <- paste0("Program ", i)
    df_plot <- merge(data@meta.data[, anno, drop = F] %>% 
                       tibble::rownames_to_column("cell"),
                     data.frame(cell = names(cur.scores),
                                score = cur.scores),
                     by = "cell") %>% 
      group_by(.dots = anno) %>% 
      summarise(mean_score = mean(score)) %>% 
      mutate(res = as.character(get({{anno}})))
    
    if (min(df_plot$mean_score, na.rm = T) < 0) {
      df_plot <- df_plot %>% arrange(mean_score) %>% 
        mutate(res = factor(res, levels = rev(unique(res)))) %>% 
        filter(mean_score != 0) %>% 
        slice_head(n = 15)
      breaks <- quantile(scales::extended_breaks()(df_plot$mean_score), probs = c(0.1, 0.5, 1))
      
      p <- ggplot(df_plot,
                  aes(x = res, y = mean_score, fill = mean_score)) +
        geom_bar(stat = "identity") +
        scale_fill_gradient(high = "grey90", low = "blue",
                            breaks = breaks, 
                            labels = breaks) +
        coord_flip() +
        scale_y_reverse() +
        theme_minimal() +
        theme(axis.title.x = element_blank(),
              axis.title.y = element_blank(),
              legend.key.size = unit(0.3, 'cm'),
              axis.text.x = element_blank(),
              axis.ticks = element_blank(),
              axis.line = element_blank(),
              panel.grid = element_blank(),
              axis.text.y = element_text(size = 15, color = "black"),
              legend.text = element_text(size = 15, angle = 0, vjust = 0, hjust = 0),
              legend.title = element_text(size = 18, vjust = 0.8),
              legend.position = "right",
              plot.title = element_text(size = 20, 
                                        hjust = 0.5)) +
        theme(plot.margin = margin(10,10,10,20)) +
        ggtitle(pr_name)
    } else {
      df_plot <- df_plot %>% arrange(-mean_score) %>% 
        mutate(res = factor(res, levels = rev(unique(res)))) %>% 
        filter(mean_score != 0) %>% 
        slice_head(n = 15)
      breaks <- quantile(scales::extended_breaks()(df_plot$mean_score), probs = c(0, 0.5, 0.9))
      p <- ggplot(df_plot,
                  aes(x = res, y = mean_score, fill = mean_score)) +
        geom_bar(stat = "identity") +
        scale_fill_gradient(low = "grey90", high = "red",
                            breaks = breaks, 
                            labels = breaks) +
        coord_flip() +
        theme_minimal() +
        theme(axis.title.x = element_blank(),
              axis.title.y = element_blank(),
              legend.key.size = unit(0.3, 'cm'),
              axis.text.x = element_blank(),
              axis.ticks = element_blank(),
              axis.line = element_blank(),
              panel.grid = element_blank(),
              axis.text.y = element_text(size = 15, color = "black"),
              legend.text = element_text(size = 15, angle = 0, vjust = 0, hjust = 0),
              legend.title = element_text(size = 18, vjust = 0.8),
              legend.position = "right",
              plot.title = element_text(size = 20, 
                                        hjust = 0.5)) +
        theme(plot.margin = margin(10,10,10,20)) +
        ggtitle(pr_name); print(p)
    }
    plots_list[[i]] <- p
  }
  return(plots_list)
}



# Modified CellChat functions

aggregateNet_my <- function(object, sources.use = NULL, targets.use = NULL, 
                            signaling = NULL, pairLR.use = NULL, remove.isolate = TRUE, 
                            thresh = 0.05, return.object = TRUE) {
  net <- object@net
  if (is.null(sources.use) & is.null(targets.use) & is.null(signaling) & is.null(pairLR.use)) {
    prob <- net$prob
    pval <- net$pval
    pval[prob == 0] <- 1
    prob[pval >= thresh] <- 0
    net$count <- apply(prob > 0, c(1,2), sum)
    net$weight <- apply(prob, c(1,2), sum)
    net$weight[is.na(net$weight)] <- 0
    net$count[is.na(net$count)] <- 0
  } else {
    df.net <- subsetCommunication(object, slot.name = "net",
                                  sources.use = sources.use, targets.use = targets.use,
                                  signaling = signaling,
                                  pairLR.use = pairLR.use,
                                  thresh = thresh)
    # View(df.net)
    df.net$source_target <- paste(df.net$source, df.net$target, sep = "|")
    df.net2 <- df.net %>% group_by(source_target) %>% summarize(count = n(), .groups = 'drop')
    df.net3 <- df.net %>% group_by(source_target) %>% summarize(prob = sum(prob), .groups = 'drop')
    df.net2$prob <- df.net3$prob
    a <- stringr::str_split(df.net2$source_target, "\\|", simplify = T)
    df.net2$source <- as.character(a[, 1, drop = F])
    df.net2$target <- as.character(a[, 2, drop = F])
    cells.level <- levels(object@idents)
    if (remove.isolate) {
      message("Isolate cell groups without any interactions are removed. To block it, set `remove.isolate = FALSE`")
      df.net2$source <- factor(df.net2$source, levels = cells.level[cells.level %in% unique(df.net2$source)])
      df.net2$target <- factor(df.net2$target, levels = cells.level[cells.level %in% unique(df.net2$target)])
    } else {
      df.net2$source <- factor(df.net2$source, levels = cells.level)
      df.net2$target <- factor(df.net2$target, levels = cells.level)
    }
    
    count <- tapply(df.net2[["count"]], list(df.net2[["source"]], df.net2[["target"]]), sum)
    prob <- tapply(df.net2[["prob"]], list(df.net2[["source"]], df.net2[["target"]]), sum)
    net$count <- count
    net$weight <- prob
    net$weight[is.na(net$weight)] <- 0
    net$count[is.na(net$count)] <- 0
  }
  if (return.object) {
    object@net <- net
    return(object)
  } else {
    return(net)
  }
  
}

renameClusters <- function(cellchat, keep_parts = 2) {
  # Extract current cluster names
  current_names <- levels(cellchat@meta$labels)
  
  # Clean names using regex to keep number and first text segment
  new_names <- sapply(current_names, function(x) {
    parts <- unlist(strsplit(x, "_")[[1]])
    paste(parts[1:min(keep_parts, length(parts))], collapse = "_")
  })
  
  # Update cluster names in all relevant slots
  levels(cellchat@idents) <- new_names
  rownames(cellchat@net[[1]]$prob) <- new_names
  colnames(cellchat@net[[1]]$prob) <- new_names
  rownames(cellchat@net[[2]]$prob) <- new_names
  colnames(cellchat@net[[2]]$prob) <- new_names
  
  rownames(cellchat@net[[1]]$count) <- new_names
  colnames(cellchat@net[[1]]$count) <- new_names
  rownames(cellchat@net[[2]]$count) <- new_names
  colnames(cellchat@net[[2]]$count) <- new_names
  
  rownames(cellchat@net[[1]]$weight) <- new_names
  colnames(cellchat@net[[1]]$weight) <- new_names
  rownames(cellchat@net[[2]]$weight) <- new_names
  colnames(cellchat@net[[2]]$weight) <- new_names
  
  # Update LR slot if needed
  if(!is.null(cellchat@LR)) {
    cellchat@LR$LRsig$cluster <- new_names[match(cellchat@LR$LRsig$cluster, current_names)]
  }
  
  return(cellchat)
}


netVisual_diffInteraction_my <- function(object, comparison = c(1,2), measure = c("count", "weight", "count.merged", "weight.merged"), color.use = NULL, color.edge = c('#b2182b','#2166ac'), title.name = NULL, sources.use = NULL, targets.use = NULL, remove.isolate = FALSE, top = 1,
                                          weight.scale = FALSE, vertex.weight = 20, vertex.weight.max = NULL, vertex.size.max = 15, vertex.label.cex=1,vertex.label.color= "black",
                                          edge.weight.max = NULL, edge.width.max=8, alpha.edge = 0.6, label.edge = FALSE,edge.label.color='black',edge.label.cex=0.8,
                                          edge.curved=0.2,shape='circle',layout=in_circle(), margin=0.2,
                                          arrow.width=1,arrow.size = 0.2){
  options(warn = -1)
  measure <- match.arg(measure)
  obj1 <- object@net[[comparison[1]]][[measure]]
  obj2 <- object@net[[comparison[2]]][[measure]]
  net.diff <- obj2 - obj1
  if (measure %in% c("count", "count.merged")) {
    if (is.null(title.name)) {
      title.name = "Differential number of interactions"
    }
  } else if (measure %in% c("weight", "weight.merged")) {
    if (is.null(title.name)) {
      title.name = "Differential interaction strength"
    }
  }
  net <- net.diff
  if ((!is.null(sources.use)) | (!is.null(targets.use))) {
    df.net <- reshape2::melt(net, value.name = "value")
    colnames(df.net)[1:2] <- c("source","target")
    # keep the interactions associated with sources and targets of interest
    if (!is.null(sources.use)){
      if (is.numeric(sources.use)) {
        sources.use <- rownames(net.diff)[sources.use]
      }
      df.net <- subset(df.net, source %in% sources.use)
    }
    if (!is.null(targets.use)){
      if (is.numeric(targets.use)) {
        targets.use <- rownames(net.diff)[targets.use]
      }
      df.net <- subset(df.net, target %in% targets.use)
    }
    cells.level <- rownames(net.diff)
    df.net$source <- factor(df.net$source, levels = cells.level)
    df.net$target <- factor(df.net$target, levels = cells.level)
    df.net$value[is.na(df.net$value)] <- 0
    net <- tapply(df.net[["value"]], list(df.net[["source"]], df.net[["target"]]), sum)
    net[is.na(net)] <- 0
  }
  
  if (remove.isolate) {
    idx1 <- which(Matrix::rowSums(net) == 0)
    idx2 <- which(Matrix::colSums(net) == 0)
    idx <- intersect(idx1, idx2)
    net <- net[-idx, ]
    net <- net[, -idx]
  }
  
  net[abs(net) < stats::quantile(abs(net), probs = 1-top, na.rm= T)] <- 0
  
  g <- graph_from_adjacency_matrix(net, mode = "directed", weighted = T)
  edge.start <- igraph::ends(g, es=igraph::E(g), names=FALSE)
  coords<-layout_(g,layout)
  if(nrow(coords)!=1){
    coords_scale=scale(coords)
  }else{
    coords_scale<-coords
  }
  if (is.null(color.use)) {
    print("Colors NULL")
    color.use = scPalette(length(igraph::V(g)))
  }
  if (is.null(vertex.weight.max)) {
    vertex.weight.max <- max(vertex.weight)
  }
  vertex.weight <- vertex.weight/vertex.weight.max*vertex.size.max+5
  
  loop.angle<-ifelse(coords_scale[igraph::V(g),1]>0,-atan(coords_scale[igraph::V(g),2]/coords_scale[igraph::V(g),1]),pi-atan(coords_scale[igraph::V(g),2]/coords_scale[igraph::V(g),1]))
  igraph::V(g)$size<-vertex.weight
  igraph::V(g)$color<-color.use[igraph::V(g)$name]
  igraph::V(g)$frame.color <- color.use[igraph::V(g)$name]
  igraph::V(g)$label.color <- vertex.label.color
  igraph::V(g)$label.cex<-vertex.label.cex
  if(label.edge){
    igraph::E(g)$label<-igraph::E(g)$weight
    igraph::E(g)$label <- round(igraph::E(g)$label, digits = 1)
  }
  igraph::E(g)$arrow.width<-arrow.width
  igraph::E(g)$arrow.size<-arrow.size
  igraph::E(g)$label.color<-edge.label.color
  igraph::E(g)$label.cex<-edge.label.cex
  #igraph::E(g)$color<- grDevices::adjustcolor(igraph::V(g)$color[edge.start[,1]],alpha.edge)
  igraph::E(g)$color <- ifelse(igraph::E(g)$weight > 0, color.edge[1],color.edge[2])
  igraph::E(g)$color <- grDevices::adjustcolor(igraph::E(g)$color, alpha.edge)
  
  igraph::E(g)$weight <- abs(igraph::E(g)$weight)
  
  if (is.null(edge.weight.max)) {
    edge.weight.max <- max(igraph::E(g)$weight)
  }
  if (weight.scale == TRUE) {
    #E(g)$width<-0.3+edge.width.max/(max(E(g)$weight)-min(E(g)$weight))*(E(g)$weight-min(E(g)$weight))
    igraph::E(g)$width<- 0.3+igraph::E(g)$weight/edge.weight.max*edge.width.max
  }else{
    igraph::E(g)$width<-0.3+edge.width.max*igraph::E(g)$weight
  }
  
  igraph::E(g)$loop.angle <- 0
  if(sum(edge.start[,2]==edge.start[,1])!=0){
    igraph::E(g)$loop.angle[which(edge.start[,2]==edge.start[,1])]<-loop.angle[edge.start[which(edge.start[,2]==edge.start[,1]),1]]
  }
  radian.rescale <- function(x, start=0, direction=1) {
    c.rotate <- function(x) (x + start) %% (2 * pi) * direction
    c.rotate(scales::rescale(x, c(0, 2 * pi), range(x)))
  }
  label.locs <- radian.rescale(x=1:length(igraph::V(g)), direction=-1, start=0)
  label.dist <- vertex.weight/max(vertex.weight)+2
  plot(g,edge.curved=edge.curved,vertex.shape=shape,layout=coords_scale,margin=margin, vertex.label.dist=label.dist,
       vertex.label.degree=label.locs, vertex.label.family="Helvetica", edge.label.family="Helvetica") # "sans"
  if (!is.null(title.name)) {
    text(0,1.5,title.name, cex = 1.1)
  }
  # https://www.andrewheiss.com/blog/2016/12/08/save-base-graphics-as-pseudo-objects-in-r/
  # grid.echo()
  # gg <-  grid.grab()
  gg <- recordPlot()
  return(gg)
}



netVisual_heatmap_my <- function(object, comparison = c(1,2), measure = c("count", "weight"), signaling = NULL, slot.name = c("netP", "net"), color.use = NULL, color.heatmap = NULL,
                                 title.name = NULL, width = NULL, height = NULL, font.size = 8, font.size.title = 10, cluster.rows = FALSE, cluster.cols = FALSE,
                                 sources.use = NULL, targets.use = NULL, remove.isolate = FALSE, row.show = NULL, col.show = NULL){
  if (!is.null(measure)) {
    measure <- match.arg(measure)
  }
  slot.name <- match.arg(slot.name)
  if (class(object@net[[1]]) == "list") {
    message("Do heatmap based on a merged object \n")
    if (is.null(color.heatmap)) {
      color.heatmap <- c('#2166ac','#b2182b')
    }
    
    # Filter for specific signaling pathway if provided
    if (!is.null(signaling)) {
      if (slot.name != "netP") {message("Set slot.name to netP")}
      prob1 <- slot(object, slot.name)[[comparison[1]]]$prob
      prob2 <- slot(object, slot.name)[[comparison[2]]]$prob
      # if (slot.name == "net") {
      #   prob[object@net$pval > thresh] <- 0
      # }
      if (!signaling %in% dimnames(prob1)[[3]]) {
        net.diff <- prob2[,,signaling]
      } else if (!signaling %in% dimnames(prob2)[[3]]) {
        net.diff <- -prob1[,,signaling]
      } else {
        net.diff <- prob2[,,signaling] - prob1[,,signaling]
      }
      
      if (is.null(title.name)) {
        title.name = paste0(signaling, " signaling network")
      }
      legend.name <- "Communication Prob. (relative)"
    } else {
      if (measure == "count") {
        if (is.null(title.name)) {
          title.name = "Differential number of interactions"
        }
      } else if (measure == "weight") {
        if (is.null(title.name)) {
          title.name = "Differential interaction strength"
        }
      }
      legend.name = "Relative values"
      
      obj1 <- object@net[[comparison[1]]][[measure]]
      obj2 <- object@net[[comparison[2]]][[measure]]
      net.diff <- obj2 - obj1
    }
    
  } else {
    message("Do heatmap based on a single object \n")
    if (is.null(color.heatmap)) {
      color.heatmap <- "Reds"
    }
    if (!is.null(signaling)) {
      if (slot.name != "netP") {message("Set slot.name to netP")}
      prob <- slot(object, slot.name)$prob
      # if (slot.name == "net") {
      #   prob[object@net$pval > thresh] <- 0
      # }
      net.diff <- prob[,,signaling]
      if (is.null(title.name)) {
        title.name = paste0(signaling, " signaling network")
      }
      legend.name <- "Communication Prob."
    } else if (!is.null(measure)) {
      net.diff <- object@net[[measure]]
      if (measure == "count") {
        if (is.null(title.name)) {
          title.name = "Number of interactions"
        }
      } else if (measure == "weight") {
        if (is.null(title.name)) {
          title.name = "Interaction strength"
        }
      }
      legend.name <- title.name
    }
  }
  
  net <- net.diff
  
  if ((!is.null(sources.use)) | (!is.null(targets.use))) {
    df.net <- reshape2::melt(net, value.name = "value")
    colnames(df.net)[1:2] <- c("source","target")
    # keep the interactions associated with sources and targets of interest
    if (!is.null(sources.use)){
      if (is.numeric(sources.use)) {
        sources.use <- rownames(net.diff)[sources.use]
      }
      df.net <- subset(df.net, source %in% sources.use)
    }
    if (!is.null(targets.use)){
      if (is.numeric(targets.use)) {
        targets.use <- rownames(net.diff)[targets.use]
      }
      df.net <- subset(df.net, target %in% targets.use)
    }
    cells.level <- rownames(net.diff)
    df.net$source <- factor(df.net$source, levels = cells.level)
    df.net$target <- factor(df.net$target, levels = cells.level)
    df.net$value[is.na(df.net$value)] <- 0
    net <- tapply(df.net[["value"]], list(df.net[["source"]], df.net[["target"]]), sum)
  }
  net[is.na(net)] <- 0
  
  if (is.null(color.use)) {
    color.use <- scPalette(ncol(net))
    names(color.use) <- colnames(net)
  }
  color.use.row <- color.use[rownames(net)]
  color.use.col <- color.use[colnames(net)]
  if (remove.isolate) {
    idx1 <- which(Matrix::rowSums(net) == 0)
    idx2 <- which(Matrix::colSums(net) == 0)
    #idx <- intersect(idx1, idx2)
    # if (length(idx) > 0) {
    #   net <- net[-idx, ]
    #   net <- net[, -idx]
    # }
    if (length(idx1) > 0) {
      net <- net[-idx1, , drop=F]
      color.use.row <- color.use.row[-idx1]
    }
    if (length(idx2) > 0) {
      net <- net[, -idx2, drop = F]
      color.use.col <- color.use.col[-idx2]
    }
    if (length(color.use.col) == 0 | length(color.use.row) == 0){
      message("After removing isolates there is no significant pairs")
      return(ggplot() + theme_minimal())
    }
    
  }
  
  mat <- net
  if (!is.null(row.show)) {
    mat <- mat[row.show, , drop=FALSE]
    color.use.row <- color.use.row[row.show]
  }
  if (!is.null(col.show)) {
    mat <- mat[ ,col.show, drop=FALSE]
    color.use.col <- color.use.col[col.show]
  }
  
  # If 1x1 matrix
  if (nrow(mat) == 1 && ncol(mat) == 1) {
    # Disable clustering and annotations for single-cell case
    cluster.rows <- FALSE
    cluster.cols <- FALSE
    ha1 <- NULL  # Remove row annotation
    ha2 <- NULL  # Remove column annotation
    
    # Ensure at least 2 distinct values for colorRamp2
    if (mat[1,1] == 0) {
      color.heatmap.use <- c("white", "#f7f7f7")  # Custom 2-color scale
    } else {
      color.heatmap.use <- circlize::colorRamp2(c(mat[1,1] - 1e-10, mat[1,1]), 
                                                c("white", ifelse(mat[1,1] < 0, '#2166ac','#b2182b')))
    }
    
    print(barplot(c(1,1), col = color.heatmap.use))
    
    # Simplify heatmap parameters
    ht1 <- Heatmap(mat, 
                   col = color.heatmap.use,
                   na_col = "white",
                   name = legend.name,
                   cluster_rows = cluster.rows,
                   cluster_columns = cluster.cols,
                   row_names_side = "left",
                   column_names_side = "bottom",
                   column_title = title.name,
                   heatmap_legend_param = list(
                     title_gp = gpar(fontsize = 8),
                     labels_gp = gpar(fontsize = 8)
                   )
    )
    return(ht1)
    
    
  } else {
    # if there positive and negative range
    if (min(mat) < 0 & max(mat > 0)) {
      color.heatmap.use = colorRamp3(c(min(mat), 0, max(mat)), c(color.heatmap[1], "#f7f7f7", color.heatmap[2]))
      colorbar.break <- c(round(min(mat, na.rm = T), digits = nchar(sub(".*\\.(0*).*","\\1",min(mat, na.rm = T)))+1), 0, round(max(mat, na.rm = T), digits = nchar(sub(".*\\.(0*).*","\\1",max(mat, na.rm = T)))+1))
    } else {
      if (length(color.heatmap) == 3) {
        color.heatmap.use = colorRamp3(c(0, min(mat), max(mat)), color.heatmap)
      } else if (length(color.heatmap) == 2) {
        # if only positive values:
        if (all(mat >= 0)) {message("All +"); color.heatmap.use = colorRamp3(c(min(mat), max(mat)), c("#f7f7f7", color.heatmap[2]))}
        # if only neative values:
        if (all(mat <= 0)) {message("All -"); color.heatmap.use = colorRamp3(c(min(mat), max(mat)), c(color.heatmap[1], "#f7f7f7"))}
        # color.heatmap.use = colorRamp3(c(min(mat), max(mat)), color.heatmap)
      } else if (length(color.heatmap) == 1) {
        color.heatmap.use = grDevices::colorRampPalette((RColorBrewer::brewer.pal(n = 9, name = color.heatmap)))(100)
      }
      colorbar.break <- c(round(min(mat, na.rm = T), digits = nchar(sub(".*\\.(0*).*","\\1",min(mat, na.rm = T)))+1), round(max(mat, na.rm = T), digits = nchar(sub(".*\\.(0*).*","\\1",max(mat, na.rm = T)))+1))
    }
    
    df.col<- data.frame(group = colnames(mat)); rownames(df.col) <- colnames(mat)
    df.row<- data.frame(group = rownames(mat)); rownames(df.row) <- rownames(mat)
    col_annotation <- HeatmapAnnotation(df = df.col, col = list(group = color.use.col),which = "column",
                                        show_legend = FALSE, show_annotation_name = FALSE,
                                        simple_anno_size = grid::unit(0.2, "cm"))
    row_annotation <- HeatmapAnnotation(df = df.row, col = list(group = color.use.row), which = "row",
                                        show_legend = FALSE, show_annotation_name = FALSE,
                                        simple_anno_size = grid::unit(0.2, "cm"))
    
    ha1 = rowAnnotation(Strength = anno_barplot(rowSums(abs(mat)), border = FALSE,gp = gpar(fill = color.use.row, col=color.use.row)), show_annotation_name = FALSE)
    ha2 = HeatmapAnnotation(Strength = anno_barplot(colSums(abs(mat)), border = FALSE,gp = gpar(fill = color.use.col, col=color.use.col)), show_annotation_name = FALSE)
    
    if (sum(abs(mat) > 0) == 1) {
      color.heatmap.use = c("white", color.heatmap.use)
    } else {
      mat[mat == 0] <- NA
    }
    
    mat[is.na(mat)] <- 0
    
    ht1 = Heatmap(mat, col = color.heatmap.use, na_col = "white", name = legend.name,
                  bottom_annotation = col_annotation, left_annotation =row_annotation, top_annotation = ha2, right_annotation = ha1,
                  cluster_rows = cluster.rows, cluster_columns = cluster.cols,
                  row_names_side = "left",row_names_rot = 0,row_names_gp = gpar(fontsize = font.size),column_names_gp = gpar(fontsize = font.size),
                  # width = unit(width, "cm"), height = unit(height, "cm"),
                  column_title = title.name,column_title_gp = gpar(fontsize = font.size.title),column_names_rot = 90,
                  row_title = "Sources (Sender)",row_title_gp = gpar(fontsize = font.size.title),row_title_rot = 90,
                  heatmap_legend_param = list(title_gp = gpar(fontsize = 8, fontface = "plain"),title_position = "leftcenter-rot",
                                              border = NA, #at = colorbar.break,
                                              legend_height = unit(20, "mm"),labels_gp = gpar(fontsize = 8),grid_width = unit(2, "mm"))
    )
    return(ht1)
  }
}


rankNet_my <- function(object, slot.name = "netP", measure = c("weight","count"), 
                       mode = c("comparison", "single"), comparison = c(1,2), 
                       color.use = NULL, stacked = FALSE, 
                       sources.use = NULL, targets.use = NULL,  
                       signaling = NULL, pairLR = NULL, signaling.type = NULL, 
                       do.stat = FALSE, paired.test = TRUE, cutoff.pvalue = 0.05, 
                       tol = 0.05, thresh = 0.05, show.raw = FALSE, return.data = FALSE, 
                       x.rotation = 90, title = NULL, bar.w = 0.75, font.size = 8,
                       do.flip = TRUE, x.angle = NULL, y.angle = 0, x.hjust = 1,y.hjust = 1,
                       axis.gap = FALSE, ylim = NULL, segments = NULL, 
                       tick_width = NULL, rel_heights = c(0.9,0,0.1)) {
  measure <- match.arg(measure)
  mode <- match.arg(mode)
  options(warn = -1)
  object.names <- names(methods::slot(object, slot.name))
  if (measure == "weight") {
    ylabel = "Information flow"
  } else if (measure == "count") {
    ylabel = "Number of interactions"
  }
  if (mode == "single") {
    object1 <- methods::slot(object, slot.name)
    prob = object1$prob
    prob[object1$pval > thresh] <- 0
    if (measure == "count") {
      prob <- 1*(prob > 0)
    }
    if (!is.null(sources.use)) {
      if (is.character(sources.use)) {
        if (all(sources.use %in% dimnames(prob)[[1]])) {
          sources.use <- match(sources.use, dimnames(prob)[[1]])
        } else {
          stop("The input `sources.use` should be cell group names or a numerical vector!")
        }
      }
      idx.t <- setdiff(1:nrow(prob), sources.use)
      prob[idx.t, , ] <- 0
    }
    if (!is.null(targets.use)) {
      if (is.character(targets.use)) {
        if (all(targets.use %in% dimnames(prob)[[1]])) {
          targets.use <- match(targets.use, dimnames(prob)[[2]])
        } else {
          stop("The input `targets.use` should be cell group names or a numerical vector!")
        }
      }
      idx.t <- setdiff(1:nrow(prob), targets.use)
      prob[ ,idx.t, ] <- 0
    }
    if (sum(prob) == 0) {
      stop("No inferred communications for the input!")
    }
    
    pSum <- apply(prob, 3, sum)
    pSum.original <- pSum
    if (measure == "weight") {
      pSum <- -1/log(pSum)
      pSum[is.na(pSum)] <- 0
      idx1 <- which(is.infinite(pSum) | pSum < 0)
      values.assign <- seq(max(pSum)*1.1, max(pSum)*1.5, length.out = length(idx1))
      position <- sort(pSum.original[idx1], index.return = TRUE)$ix
      pSum[idx1] <- values.assign[match(1:length(idx1), position)]
    } else if (measure == "count") {
      pSum <- pSum.original
    }
    
    pair.name <- names(pSum)
    
    df<- data.frame(name = pair.name, contribution = pSum.original, contribution.scaled = pSum, group = object.names[comparison[1]])
    idx <- with(df, order(df$contribution))
    df <- df[idx, ]
    df$name <- factor(df$name, levels = as.character(df$name))
    print(dim(df)) # !!!!
    for (i in 1:length(pair.name)) {
      df.t <- df[df$name == pair.name[i], "contribution"]
      if (sum(df.t) == 0) {
        df <- df[-which(df$name == pair.name[i]), ]
      }
    }
    print(dim(df)) # !!!!
    
    if (!is.null(signaling.type)) {
      LR <- subset(object@DB$interaction, annotation %in% signaling.type)
      if (slot.name == "netP") {
        signaling <- unique(LR$pathway_name)
      } else if (slot.name == "net") {
        pairLR <- LR$interaction_name
      }
    }
    
    if ((slot.name == "netP") && (!is.null(signaling))) {
      df <- subset(df, name %in% signaling)
    } else if ((slot.name == "netP") &&(!is.null(pairLR))) {
      stop("You need to set `slot.name == 'net'` if showing specific L-R pairs ")
    }
    if ((slot.name == "net") && (!is.null(pairLR))) {
      df <- subset(df, name %in% pairLR)
      print(dim(df)) # !!!!
    } else if ((slot.name == "net") && (!is.null(signaling))) {
      stop("You need to set `slot.name == 'netP'` if showing specific signaling pathways ")
    }
    
    gg <- ggplot(df, aes(x=name, y=contribution.scaled)) + geom_bar(stat="identity",width = bar.w) +
      theme_classic() + theme(axis.text=element_text(size=font.size),axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.y = element_text(size=10)) +
      xlab("") + ylab(ylabel) + coord_flip()#+
    if (!is.null(title)) {
      gg <- gg + ggtitle(title)+ theme(plot.title = element_text(hjust = 0.5))
    }
    
  } else if (mode == "comparison") {
    prob.list <- list()
    pSum <- list()
    pSum.original <- list()
    pair.name <- list()
    idx <- list()
    pSum.original.all <- c()
    object.names.comparison <- c()
    for (i in 1:length(comparison)) {
      object.list <- methods::slot(object, slot.name)[[comparison[i]]]
      prob <- object.list$prob
      prob[object.list$pval > thresh] <- 0
      if (measure == "count") {
        prob <- 1*(prob > 0)
      }
      if (!is.null(sources.use)) {
        if (is.character(sources.use)) {
          if (all(sources.use %in% dimnames(prob)[[1]])) {
            sources.use <- match(sources.use, dimnames(prob)[[1]])
          } else {
            stop("The input `sources.use` should be cell group names or a numerical vector!")
          }
        }
        idx.t <- setdiff(1:nrow(prob), sources.use)
        prob[idx.t, , ] <- 0
      }
      if (!is.null(targets.use)) {
        if (is.character(targets.use)) {
          if (all(targets.use %in% dimnames(prob)[[1]])) {
            targets.use <- match(targets.use, dimnames(prob)[[2]])
          } else {
            stop("The input `targets.use` should be cell group names or a numerical vector!")
          }
        }
        idx.t <- setdiff(1:nrow(prob), targets.use)
        prob[ ,idx.t, ] <- 0
      }
      if (sum(prob) == 0) {
        stop("No inferred communications for the input!")
      }
      prob.list[[i]] <- prob
      pSum.original[[i]] <- apply(prob, 3, sum)
      if (measure == "weight") {
        pSum[[i]] <- -1/log(pSum.original[[i]])
        pSum[[i]][is.na(pSum[[i]])] <- 0
        idx[[i]] <- which(is.infinite(pSum[[i]]) | pSum[[i]] < 0)
        pSum.original.all <- c(pSum.original.all, pSum.original[[i]][idx[[i]]])
      } else if (measure == "count") {
        pSum[[i]] <- pSum.original[[i]] # the prob is already binarized in line 1136
      }
      pair.name[[i]] <- names(pSum.original[[i]])
      object.names.comparison <- c(object.names.comparison, object.names[comparison[i]])
    }
    if (measure == "weight") {
      values.assign <- seq(max(unlist(pSum))*1.1, max(unlist(pSum))*1.5, length.out = length(unlist(idx)))
      position <- sort(pSum.original.all, index.return = TRUE)$ix
      for (i in 1:length(comparison)) {
        if (i == 1) {
          pSum[[i]][idx[[i]]] <- values.assign[match(1:length(idx[[i]]), position)]
        } else {
          pSum[[i]][idx[[i]]] <- values.assign[match(length(unlist(idx[1:i-1]))+1:length(unlist(idx[1:i])), position)]
        }
      }
    }
    
    
    
    pair.name.all <- as.character(unique(unlist(pair.name)))
    df <- list()
    for (i in 1:length(comparison)) {
      df[[i]] <- data.frame(name = pair.name.all, contribution = 0, contribution.scaled = 0, group = object.names[comparison[i]], row.names = pair.name.all)
      df[[i]][pair.name[[i]],3] <- pSum[[i]]
      df[[i]][pair.name[[i]],2] <- pSum.original[[i]]
    }
    
    
    # contribution.relative <- as.numeric(format(df[[length(comparison)]]$contribution/abs(df[[1]]$contribution), digits=1))
    # #  contribution.relative <- as.numeric(format(df[[length(comparison)]]$contribution.scaled/abs(df[[1]]$contribution.scaled), digits=1))
    # contribution.relative2 <- as.numeric(format(df[[length(comparison)-1]]$contribution/abs(df[[1]]$contribution), digits=1))
    # contribution.relative[is.na(contribution.relative)] <- 0
    # for (i in 1:length(comparison)) {
    #   df[[i]]$contribution.relative <- contribution.relative
    #   df[[i]]$contribution.relative2 <- contribution.relative2
    # }
    # df[[1]]$contribution.data2 <- df[[length(comparison)]]$contribution
    # idx <- with(df[[1]], order(-contribution.relative,  -contribution.relative2, contribution, -contribution.data2))
    #
    contribution.relative <- list()
    for (i in 1:(length(comparison)-1)) {
      contribution.relative[[i]] <- as.numeric(format(df[[length(comparison)-i+1]]$contribution/df[[1]]$contribution, digits=1))
      contribution.relative[[i]][is.na(contribution.relative[[i]])] <- 0
    }
    names(contribution.relative) <- paste0("contribution.relative.", 1:length(contribution.relative))
    for (i in 1:length(comparison)) {
      for (j in 1:length(contribution.relative)) {
        df[[i]][[names(contribution.relative)[j]]] <- contribution.relative[[j]]
      }
    }
    df[[1]]$contribution.data2 <- df[[length(comparison)]]$contribution
    if (length(comparison) == 2) {
      idx <- with(df[[1]], order(-contribution.relative.1, contribution, -contribution.data2))
    } else if (length(comparison) == 3) {
      idx <- with(df[[1]], order(-contribution.relative.1, -contribution.relative.2,contribution, -contribution.data2))
    } else if (length(comparison) == 4) {
      idx <- with(df[[1]], order(-contribution.relative.1, -contribution.relative.2, -contribution.relative.3, contribution, -contribution.data2))
    } else {
      idx <- with(df[[1]], order(-contribution.relative.1, -contribution.relative.2, -contribution.relative.3, -contribution.relative.4, contribution, -contribution.data2))
    }
    
    
    
    for (i in 1:length(comparison)) {
      df[[i]] <- df[[i]][idx, ]
      df[[i]]$name <- factor(df[[i]]$name, levels = as.character(df[[i]]$name))
    }
    df[[1]]$contribution.data2 <- NULL
    
    df <- do.call(rbind, df)
    df$group <- factor(df$group, levels = object.names.comparison)
    
    if (is.null(color.use)) {
      color.use =  ggPalette(length(comparison))
    }
    
    # https://stackoverflow.com/questions/49448497/coord-flip-changes-ordering-of-bars-within-groups-in-grouped-bar-plot
    df$group <- factor(df$group, levels = rev(levels(df$group)))
    color.use <- rev(color.use)
    
    # perform statistical analysis
    # if (do.stat) {
    #   pvalues <- c()
    #   for (i in 1:length(pair.name.all)) {
    #     df.prob <- data.frame()
    #     for (j in 1:length(comparison)) {
    #       if (pair.name.all[i] %in% pair.name[[j]]) {
    #         df.prob <- rbind(df.prob, data.frame(prob = as.vector(prob.list[[j]][ , , pair.name.all[i]]), group = comparison[j]))
    #       } else {
    #         df.prob <- rbind(df.prob, data.frame(prob = as.vector(matrix(0, nrow = nrow(prob.list[[j]]), ncol = nrow(prob.list[[j]]))), group = comparison[j]))
    #       }
    #
    #     }
    #     df.prob$group <- factor(df.prob$group, levels = comparison)
    #     if (length(comparison) == 2) {
    #       pvalues[i] <- wilcox.test(prob ~ group, data = df.prob)$p.value
    #     } else {
    #       pvalues[i] <- kruskal.test(prob ~ group, data = df.prob)$p.value
    #     }
    #   }
    #   df$pvalues <- pvalues
    # }
    if (do.stat & length(comparison) == 2) {
      for (i in 1:length(pair.name.all)) {
        if (nrow(prob.list[[j]]) != nrow(prob.list[[1]])) {
          if (paired.test) {
            stop("Paired test is not applicable to datasets with different cellular compositions! Please set `do.stat = FALSE` or `paired.test = FALSE`! \n")
          }
        }
        prob.values <- matrix(0, nrow = nrow(prob.list[[1]]) * nrow(prob.list[[1]]), ncol = length(comparison))
        for (j in 1:length(comparison)) {
          if (pair.name.all[i] %in% pair.name[[j]]) {
            prob.values[, j] <- as.vector(prob.list[[j]][ , , pair.name.all[i]])
          } else {
            prob.values[, j] <- NA
          }
        }
        prob.values <- prob.values[rowSums(prob.values, na.rm = TRUE) != 0, , drop = FALSE]
        if (nrow(prob.values) >3 & sum(is.na(prob.values)) == 0) {
          pvalues <- wilcox.test(prob.values[ ,1], prob.values[ ,2], paired = paired.test)$p.value
        } else {
          pvalues <- 0
        }
        pvalues[is.na(pvalues)] <- 0
        df$pvalues[df$name == pair.name.all[i]] <- pvalues
      }
    }
    
    
    # if (length(comparison) == 2) {
    #   print(dim(df))
    #   if (do.stat) {
    #     colors.text <- ifelse((df$contribution.relative < 1-tol) & (df$pvalues < cutoff.pvalue), color.use[2], ifelse((df$contribution.relative > 1+tol) & df$pvalues < cutoff.pvalue, color.use[1], "black"))
    #   } else {
    #     colors.text <- ifelse(df$contribution.relative < 1-tol, color.use[2], ifelse(df$contribution.relative > 1+tol, color.use[1], "black"))
    #   }
    # } else {
    #   message("The text on the y-axis will not be colored for the number of compared datasets larger than 3!")
    #   colors.text = NULL
    # }
    
    for (i in 1:length(pair.name.all)) {
      df.t <- df[df$name == pair.name.all[i], "contribution"]
      if (sum(df.t) == 0) {
        df <- df[-which(df$name == pair.name.all[i]), ]
      }
    }
    
    if ((slot.name == "netP") && (!is.null(signaling))) {
      df <- subset(df, name %in% signaling)
    } else if ((slot.name == "netP") &&(!is.null(pairLR))) {
      stop("You need to set `slot.name == 'net'` if showing specific L-R pairs ")
    }
    if ((slot.name == "net") && (!is.null(pairLR))) {
      df <- subset(df, name %in% pairLR)
    } else if ((slot.name == "net") && (!is.null(signaling))) {
      stop("You need to set `slot.name == 'netP'` if showing specific signaling pathways ")
    }
    
    if (length(comparison) == 2) {
      # print(dim(df))
      if (do.stat) {
        colors.text <- ifelse((df$contribution.relative < 1-tol) & (df$pvalues < cutoff.pvalue), color.use[2], ifelse((df$contribution.relative > 1+tol) & df$pvalues < cutoff.pvalue, color.use[1], "black"))
      } else {
        colors.text <- ifelse(df$contribution.relative < 1-tol, color.use[2], ifelse(df$contribution.relative > 1+tol, color.use[1], "black"))
      }
    } else {
      message("The text on the y-axis will not be colored for the number of compared datasets larger than 3!")
      colors.text = NULL
    }
    
    if (stacked) {
      gg <- ggplot(df, aes(x=name, y=contribution, fill = group)) + geom_bar(stat="identity",width = bar.w, position ="fill") # +
      # xlab("") + ylab("Relative information flow") #+ theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())
      #  scale_y_discrete(breaks=c("0","0.5","1")) +
      if (measure == "weight") {
        gg <- gg + xlab("") + ylab("Relative information flow")
      } else if (measure == "count") {
        gg <- gg + xlab("") + ylab("Relative number of interactions")
      }
      
      gg <- gg + geom_hline(yintercept = 0.5, linetype="dashed", color = "grey50", size=0.5)
    } else {
      if (show.raw) {
        gg <- ggplot(df, aes(x=name, y=contribution, fill = group)) + geom_bar(stat="identity",width = bar.w, position = position_dodge(0.8)) +
          xlab("") + ylab(ylabel) #+ coord_flip()#+ theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())
      } else {
        gg <- ggplot(df, aes(x=name, y=contribution.scaled, fill = group)) + geom_bar(stat="identity",width = bar.w, position = position_dodge(0.8)) +
          xlab("") + ylab(ylabel) #+ coord_flip()#+ theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())
      }
      
      if (axis.gap) {
        gg <- gg + theme_bw() + theme(panel.grid = element_blank())
        gg.gap::gg.gap(gg,
                       ylim = ylim,
                       segments = segments,
                       tick_width = tick_width,
                       rel_heights = rel_heights)
      }
    }
    gg <- gg +  CellChat_theme_opts() + theme_classic()
    if (do.flip) {
      gg <- gg + coord_flip() + theme(axis.text.y = element_text(colour = colors.text))
      if (is.null(x.angle)) {
        x.angle = 0
      }
      
    } else {
      if (is.null(x.angle)) {
        x.angle = 45
      }
      # print(colors.text)
      gg <- gg + scale_x_discrete(limits = rev) + theme(axis.text.x = element_text(colour = rev(colors.text)))
      
    }
    
    gg <- gg + theme(axis.text=element_text(size=font.size), axis.title.y = element_text(size=font.size))
    gg <- gg + scale_fill_manual(name = "", values = color.use)
    gg <- gg + guides(fill = guide_legend(reverse = TRUE))
    gg <- gg + theme(axis.text.x = element_text(angle = x.angle, hjust=x.hjust),
                     axis.text.y = element_text(angle = y.angle, hjust=y.hjust))
    if (!is.null(title)) {
      gg <- gg + ggtitle(title)+ theme(plot.title = element_text(hjust = 0.5))
    }
  }
  
  if (return.data) {
    df$contribution <- abs(df$contribution)
    df$contribution.scaled <- abs(df$contribution.scaled)
    return(list(signaling.contribution = df, gg.obj = gg))
  } else {
    return(gg)
  }
}



volcano_de_genes <- function(markers, 
                             p_val_adj_cutoff=0.05, 
                             avg_log2FC_cutoff=0.25,
                             pct_cutoff=0.1,
                             genes_highlight = NULL,
                             txt.size = 4,
                             label.left = F,
                             up_label, down_label,
                             ...) {

  markers$p_val_adj <- ifelse(markers$p_val_adj == 0, 
                              1e-323,
                              markers$p_val_adj)
  markers$pct.diff <- markers$pct.1 - markers$pct.2
  
  markers2 <- markers %>% filter(p_val_adj < p_val_adj_cutoff, 
                                 abs(avg_log2FC) > avg_log2FC_cutoff,
                                 pct.1 > pct_cutoff | pct.2 > pct_cutoff)
  
  
  up <- (markers2 %>% filter(avg_log2FC > 0, gene %in% genes_highlight))$gene
  down <- (markers2 %>% filter(avg_log2FC < 0, gene %in% genes_highlight))$gene
  
  top <- union(up, down)
  
  markers_plot <- markers %>%
    mutate(diff_expr = case_when(
      gene %in% up ~ up_label,
      gene %in% down ~ down_label,
      TRUE ~ "other"
    ),
    label = ifelse(gene %in% top, gene, NA))
  
  markers_plot$diff_expr <- factor(markers_plot$diff_expr,
                                   levels = rev(c(up_label, down_label, "other")))
  
  colors_diff_expr <- c("red", "blue", "grey55")
  names(colors_diff_expr) <- c(up_label, down_label, "other")
  colors_diff_expr <- colors_diff_expr[names(colors_diff_expr) %in% unique(markers_plot$diff_expr)]
  
  if (label.left) {
    ggplot(data=markers_plot %>% 
             arrange(diff_expr), 
           aes(x=avg_log2FC, y=-log10(p_val_adj), label=label)) +
      ggnewscale::new_scale_color() +
      geom_point(aes(alpha = abs(pct.diff), col=diff_expr)) +
      scale_color_manual(values=colors_diff_expr) +
      theme_minimal() +
      ggrepel::geom_text_repel(aes(color = diff_expr), max.overlaps = 100, 
                               size = txt.size, show.legend = FALSE, 
                               xlim = c(min(markers_plot$avg_log2FC), 0),  ...) +
      geom_vline(xintercept=c(-avg_log2FC_cutoff, avg_log2FC_cutoff), 
                 col="grey30", linetype = "dashed") +
      geom_hline(yintercept=-log10(p_val_adj_cutoff), col="grey30", linetype = "dashed") +
      theme(plot.margin = unit(c(1, 0, 0, 0), "cm")) +
      guides(alpha=guide_legend(title="pct.diff"))
  } else {
    ggplot(data=markers_plot %>% 
             arrange(diff_expr), 
           aes(x=avg_log2FC, y=-log10(p_val_adj), label=label)) +
      ggnewscale::new_scale_color() +
      geom_point(aes(alpha = abs(pct.diff), col=diff_expr)) +
      scale_color_manual(values=colors_diff_expr) +
      theme_minimal() +
      ggrepel::geom_text_repel(aes(color = diff_expr), max.overlaps = 100, 
                               size = txt.size, show.legend = FALSE) +
      geom_vline(xintercept=c(-avg_log2FC_cutoff, avg_log2FC_cutoff), 
                 col="grey30", linetype = "dashed") +
      geom_hline(yintercept=-log10(p_val_adj_cutoff), col="grey30", linetype = "dashed") +
      theme(plot.margin = unit(c(1, 0, 0, 0), "cm")) +
      guides(alpha=guide_legend(title="pct.diff"))
  }
  
  
}



swap_elements <- function(v, x, y) {
  vx <- v[x]
  v[x] <- v[y]
  v[y] <- vx
  v
}