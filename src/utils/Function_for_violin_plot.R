# ==============================================================================
# Helper functions for violin plots and statistical summaries
# ==============================================================================
# Source: https://gist.github.com/Karel-Kroeze/746685f5613e01ba820a31e57f87ec87
# ==============================================================================

# ------------------------------------------------------------------------------
# summarySE: Calculate summary statistics with error bars
# ------------------------------------------------------------------------------
#' @param data Data frame containing the data
#' @param measurevar Name of the column containing the measurement variable
#' @param groupvars Vector of column names to group by
#' @param na.rm Boolean to remove NA values
#' @param conf.interval Confidence interval (default 0.95)
#' @param .drop Whether to drop unused factor levels
#' @return Data frame with N, mean, sd, se, and ci columns
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)

  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }

  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )

  # Rename the "mean" column
  datac <- plyr::rename(datac, c("mean" = measurevar))

  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean

  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval:
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult

  return(datac)
}

# ------------------------------------------------------------------------------
# GeomSplitViolin: Custom ggproto for split violin plots
# ------------------------------------------------------------------------------
GeomSplitViolin <- ggproto(
  "GeomSplitViolin",
  GeomViolin,
  draw_group = function(self, data, ..., draw_quantiles = NULL) {
    data <- base::transform(data,
                            xminv = x - violinwidth * (x - xmin),
                            xmaxv = x + violinwidth * (xmax - x))
    grp <- data[1,'group']
    newdata <- plyr::arrange(
      base:: transform(data, x = if(grp%%2==1) xminv else xmaxv),
      if(grp%%2==1) y else -y
    )
    newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
    newdata[c(1,nrow(newdata)-1,nrow(newdata)), 'x'] <- round(newdata[1, 'x'])
    if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
      stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <= 1))
      quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
      aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
      aesthetics$alpha <- rep(1, nrow(quantiles))
      both <- cbind(quantiles, aesthetics)
      quantile_grob <- GeomPath$draw_panel(both, ...)
      ggplot2:::ggname("geom_split_violin",
                       grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
    } else {
      ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
    }
  }
)

# ------------------------------------------------------------------------------
# geom_split_violin: Layer function for split violin plots
# ------------------------------------------------------------------------------
#' @param mapping Aesthetic mapping
#' @param data Data frame
#' @param stat Statistical transformation
#' @param position Position adjustment
#' @param draw_quantiles Quantiles to draw
#' @param trim Whether to trim the violin to the data range
#' @param scale Scaling method for violin widths
#' @param na.rm Remove NA values
#' @param show.legend Show legend
#' @param inherit.aes Inherit aesthetics from plot
geom_split_violin <- function (mapping = NULL,
                               data = NULL,
                               stat = "ydensity",
                               position = "identity", ...,
                               draw_quantiles = NULL,
                               trim = TRUE,
                               scale = "area",
                               na.rm = FALSE,
                               show.legend = NA,
                               inherit.aes = TRUE) {
  layer(data = data,
        mapping = mapping,
        stat = stat,
        geom = GeomSplitViolin,
        position = position,
        show.legend = show.legend,
        inherit.aes = inherit.aes,
        params = list(trim = trim,
                      scale = scale,
                      draw_quantiles = draw_quantiles,
                      na.rm = na.rm, ...)
  )
}
