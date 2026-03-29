# modified version from simplifyEnrichment
# Secondary heatmap colorscale modified to blue shades from green and red scale
simplifyGOFromMultipleLists <- function (lt, go_id_column = NULL, padj_column = NULL, padj_cutoff = 0.01, 
          filter = function(x) any(x < padj_cutoff), default = 1, ont = NULL, 
          db = "org.Hs.eg.db", measure = "Rel", heatmap_param = list(NULL), 
          show_barplot = TRUE, method = "binary_cut", control = list(), 
          min_term = NULL, verbose = TRUE, column_title = NULL, ...) 
{
  n = length(lt)
  if (is.data.frame(lt[[1]])) {
    if (is.null(go_id_column)) {
      go_id_column = which(sapply(lt[[1]], function(x) all(grepl("^GO:\\d+$", 
                                                                 x))))[1]
      if (length(go_id_column) == 0) {
        if (!is.null(rownames(lt[[1]]))) {
          go_id_column = rownames
          if (is.null(rownames(lt[[1]]))) {
            stop_wrap("Cannot find the GO ID column in the data frames. Please explicitly set argument `go_id_column`.")
          }
          if (verbose) {
            simplifyEnrichment:::message_wrap("Use row names of the data frame as `go_id_column`.")
          }
        }
        else {
          stop_wrap("Cannot find the GO ID column in the data frames. Please explicitly set argument `go_id_column`.")
        }
      }
      else {
        if (verbose) {
          simplifyEnrichment:::message_wrap(GetoptLong::qq("Use column '@{colnames(lt[[1]])[go_id_column]}' as `go_id_column`."))
        }
      }
    }
    if (is.null(padj_column)) {
      cn = colnames(lt[[1]])
      ind = simplifyEnrichment:::test_padj_column(cn)
      if (length(ind)) {
        padj_column = ind
        if (verbose) {
          simplifyEnrichment:::message_wrap(GetoptLong::qq("Use column '@{colnames(lt[[1]])[padj_column]}' as `padj_column`."))
        }
      }
      else {
        stop_wrap("Cannot find the column the contains adjusted p-values in the data frames. Please explicitly set argument `padj_column`.")
      }
    }
    lt = lapply(lt, function(x) {
      if (is.function(go_id_column)) {
        structure(x[, padj_column], names = go_id_column(x))
      }
      else {
        structure(x[, padj_column], names = x[, go_id_column])
      }
    })
    return(simplifyGOFromMultipleLists(lt, padj_cutoff = padj_cutoff, 
                                       filter = filter, default = default, ont = ont, db = db, 
                                       measure = measure, heatmap_param = heatmap_param, 
                                       show_barplot = show_barplot, method = method, control = control, 
                                       min_term = min_term, verbose = verbose, column_title = column_title, 
                                       ...))
  }
  else if (is.character(lt[[1]])) {
    lt = lapply(lt, function(x) structure(rep(1, length(x)), 
                                          names = x))
    return(simplifyGOFromMultipleLists(lt, default = 0, filter = function(x) TRUE, 
                                       ont = ont, db = db, measure = measure, show_barplot = show_barplot, 
                                       method = method, heatmap_param = list(transform = function(x) x, 
                                                                             breaks = c(0, 1), col = c("transparent", "red"), 
                                                                             name = "", labels = c("not available", "available")), 
                                       control = control, min_term = min_term, verbose = verbose, 
                                       column_title = column_title, ...))
  }
  heatmap_param2 = list(transform = NULL, breaks = NULL, col = NULL, 
                        labels = NULL, name = "padj")
  for (nm in names(heatmap_param)) {
    heatmap_param2[[nm]] = heatmap_param[[nm]]
  }
  transform = heatmap_param2$transform
  if (is.null(transform)) 
    transform = function(x) -log10(x)
  breaks = heatmap_param2$breaks
  col = heatmap_param2$col
  labels = heatmap_param2$labels
  name = heatmap_param2$name
  if (is.null(name)) 
    name = ""
  if (is.null(breaks) && is.null(col)) {
    digit = ceiling(-log10(padj_cutoff))
    base = padj_cutoff * 10^digit
    breaks = c(1, padj_cutoff, base * 10^(-digit * 2))
    col = c('white', 'cyan', 'blue4')
    labels = gt_render(c("1", GetoptLong::qq("@{base}x10<sup>-@{digit}</sup>"), 
                         GetoptLong::qq("@{base}x10<sup>-@{digit*2}</sup>")))
  }
  else if (!is.null(breaks) && !is.null(col)) {
    if (length(breaks) != length(col)) {
      stop_wrap("Length of `breaks` must be the same as the length of `col`.")
    }
  }
  all_go_id = unique(unlist(lapply(lt, names)))
  if (!all(grepl("^GO:\\d+$", all_go_id))) {
    stop_wrap("Only GO ID is allowed.")
  }
  m = matrix(default, nrow = length(all_go_id), ncol = n)
  rownames(m) = all_go_id
  colnames(m) = names(lt)
  if (is.null(colnames)) 
    colnames = paste0("Group", 1:n)
  for (i in 1:n) {
    m[names(lt[[i]]), i] = lt[[i]]
  }
  l = apply(m, 1, function(x) {
    if (all(is.na(x))) {
      FALSE
    }
    else {
      l = filter(x[!is.na(x)])
      if (length(l) == 1) {
        return(l)
      }
      else {
        return(any(l))
      }
    }
  })
  m = m[l, , drop = FALSE]
  m = t(apply(m, 1, transform))
  if (verbose) 
    message(GetoptLong::qq("@{nrow(m)}/@{length(all_go_id)} GO IDs left for clustering."))
  if (length(unique(m[!is.na(m)])) <= 2) {
    col = structure(col, names = breaks)
  }
  else {
    if (is.null(breaks) && is.null(col)) {
      col = NULL
    }
    else if (!is.null(breaks) && !is.null(col)) {
      if (length(breaks) != length(col)) {
        stop_wrap("Length of `breaks` and `col` should be the same.")
      }
      col = colorRamp2(transform(breaks), col)
    }
    else {
      stop_wrap("Arguments `breaks` and `col` should be set at the same time.")
    }
  }
  all_go_id = rownames(m)
  sim_mat = GO_similarity(all_go_id, ont = ont, db = db, measure = measure)
  all_go_id = rownames(sim_mat)
  heatmap_legend_param = list()
  heatmap_legend_param$at = transform(breaks)
  heatmap_legend_param$labels = if (is.null(labels)) 
    breaks
  else labels
  heatmap_legend_param$title = name
  mm = m[all_go_id, , drop = FALSE]
  if (show_barplot) {
    draw_ht = function(align_to) {
      s = sapply(align_to, function(index) max(apply(mm[index, 
      ], 2, function(x) sum(x >= transform(padj_cutoff)))))
      max = max(s)
      by = diff(grid.pretty(c(0, max)))[1]
      Heatmap(mm, col = col, name = if (name == "") 
        NULL
        else name, show_row_names = FALSE, cluster_columns = FALSE, 
        border = "black", heatmap_legend_param = heatmap_legend_param, 
        width = unit(0.5, "cm") * n, use_raster = TRUE, 
        left_annotation = rowAnnotation(empty = anno_block(width = unit(1.2, 
                                                                        "cm"), panel_fun = function(index) grid.text(GetoptLong::qq("Number of significant GO terms in each cluster (padj < @{padj_cutoff})"), 
                                                                                                                     unit(0, "npc"), 0.5, just = "top", rot = 90, 
                                                                                                                     gp = gpar(fontsize = 10))), bar = anno_link(align_to = align_to, 
                                                                                                                                                                 side = "left", gap = unit(3, "mm"), link_gp = gpar(fill = "#DDDDDD", 
                                                                                                                                                                                                                    col = "#AAAAAA"), internal_line = FALSE, 
                                                                                                                                                                 panel_fun = function(index) {
                                                                                                                                                                   v = apply(mm[index, ], 2, function(x) sum(x >= 
                                                                                                                                                                                                               transform(padj_cutoff)))
                                                                                                                                                                   grid.text(v[2])
                                                                                                                                                                   pushViewport(viewport())
                                                                                                                                                                   grid.rect(gp = gpar(fill = "#DDDDDD", col = "#DDDDDD"))
                                                                                                                                                                   grid.lines(c(1, 0, 0, 1), c(0, 0, 1, 1), 
                                                                                                                                                                              gp = gpar(col = "#AAAAAA"), default.units = "npc")
                                                                                                                                                                   pushViewport(viewport(xscale = c(0.5, length(v) + 
                                                                                                                                                                                                      0.5), yscale = c(0, max(v)), height = unit(1, 
                                                                                                                                                                                                                                                 "npc") - unit(2, "mm")))
                                                                                                                                                                   grid.rect(seq_along(v), 0, width = 0.6, height = unit(v, 
                                                                                                                                                                                                                         "native"), default.units = "native", just = "bottom", 
                                                                                                                                                                             gp = gpar(fill = "#444444", col = "#444444"))
                                                                                                                                                                   if (length(index)/nrow(mm) > 0.05) {
                                                                                                                                                                     grid.yaxis(at = seq(0, max(v), by = by), 
                                                                                                                                                                                gp = gpar(col = "#444444", cex = 0.6))
                                                                                                                                                                   }
                                                                                                                                                                   popViewport()
                                                                                                                                                                   popViewport()
                                                                                                                                                                 }, size = s/sum(s) * (unit(1, "npc") - unit(3, 
                                                                                                                                                                                                             "mm") * (length(align_to) - 1) - unit(2, 
                                                                                                                                                                                                                                                   "mm") * length(align_to)) + unit(2, "mm"))), 
        post_fun = function(ht) {
          decorate_annotation("bar", {
            nc = ncol(mm)
            grid.text(colnames(mm), (seq_len(nc) - 0.5)/nc * 
                        (unit(1, "npc") - unit(5, "mm")), y = -ht_opt$COLUMN_ANNO_PADDING, 
                      default.units = "npc", just = "right", 
                      rot = 90)
          })
        })
    }
  }
  else {
    draw_ht = Heatmap(mm, col = col, name = if (name == "") 
      NULL
      else name, show_row_names = FALSE, cluster_columns = FALSE, 
      border = "black", heatmap_legend_param = heatmap_legend_param, 
      width = unit(0.5, "cm") * n, use_raster = TRUE)
  }
  if (is.null(min_term)) 
    min_term = round(nrow(sim_mat) * 0.02)
  if (is.null(column_title)) 
    column_title = GetoptLong::qq("@{length(all_go_id)} GO terms clustered by '@{method}'")
  simplifyGO(sim_mat, ht_list = draw_ht, method = method, verbose = verbose, 
             min_term = min_term, control = control, column_title = column_title, 
             ...)
}