fviz_pca_biplot2 <- function (X, axes = c(1, 2), geom = c("point", "text"), geom.ind = geom, 
          geom.var = c("arrow", "text"), col.ind = "black", fill.ind = "white", 
          col.var = "steelblue", fill.var = "white", gradient.cols = NULL, 
          label = "all", invisible = "none", repel = FALSE, habillage = "none", 
          palette = NULL, addEllipses = FALSE, title = "PCA - Biplot", 
          ...) 
{
  is.individuals.colored.by.variable <- FALSE
  is.variables.colored.by.variable <- FALSE
  is.gradient.color <- FALSE
  is.gradient.fill <-FALSE
  is.discrete.color <- FALSE
  is.discrete.fill <- FALSE
  var <- facto_summarize(X, element = "var", result = c("coord", 
                                                        "contrib", "cos2"), axes = axes)
  colnames(var)[2:3] <- c("x", "y")
  pca.ind <- get_pca_ind(X)
  ind <- data.frame(pca.ind$coord[, axes, drop = FALSE])
  colnames(ind) <- c("x", "y")
  r <- min((max(ind[, "x"]) - min(ind[, "x"])/(max(var[, "x"]) - 
                                                 min(var[, "x"]))), (max(ind[, "y"]) - min(ind[, "y"])/(max(var[, 
                                                                                                                "y"]) - min(var[, "y"]))))
  ellipse.border.remove <- FALSE
  if (is.individuals.colored.by.variable & is.variables.colored.by.variable) 
    ellipse.border.remove <- TRUE
  p <- fviz_pca_ind(X, axes = axes, geom = geom.ind, repel = repel, 
                    col.ind = col.ind, fill.ind = fill.ind, label = label, 
                    invisible = invisible, habillage = habillage, addEllipses = addEllipses, 
                    ellipse.border.remove = ellipse.border.remove, ...)
  p <- fviz_pca_var(X, axes = axes, geom = geom.var, repel = repel, 
                    col.var = col.var, fill.var = fill.var, label = label, 
                    invisible = invisible, scale. = r * 0.7, ggp = p, ...)
  if (!is.null(gradient.cols)) {
    if (is.gradient.color) 
      p <- p + ggpubr::gradient_color(gradient.cols)
    if (is.gradient.fill) 
      p <- p + ggpubr::gradient_fill(gradient.cols)
  }
  if (!is.null(palette)) {
    if (is.discrete.color) 
      p <- p + ggpubr::color_palette(palette)
    if (is.discrete.fill) 
      p <- p + ggpubr::fill_palette(palette)
  }
  p + labs(title = title)
}
