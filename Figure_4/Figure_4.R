# Script for Figure 4, applying the code from
# Margvelani et al., Micro-RNAs, their target proteins, predispositions and
# the memory of filial imprinting. Scientific Reports volume 8, Article number:
# 17444 (2018)
# https://www.nature.com/articles/s41598-018-35097-w
#
# All genes, sides and regions are analysed; results are written in
# Results/<gene>/<side>_<region>/ (figure.png and stats.csv).
# Panels of Figure 4:
# a) RORA, Left IMM
# b) ROBO1, Left IMM
# c) LUC7L, Left IMM
# d) FOXP2, Left IMM
# e) FOXP2, Left PN
# f) GLUBK89, Left IMM
# g) ENSGALG00010026609, Left IMM

# set up
rm(list = ls())
source('../ancillary/StatPlot.R')
source('../ancillary/figure_settings.R')
library(readxl)

# control panel
data_folder <- '../data/figure_4'
sides <- c('L', 'R')
regions <- c('IMM', 'PN')
folders <- dir(data_folder)
root_res_folder <- 'Results'

# panels of Figure 4: (a | b) / (c | d | e) / (f | g), sizes in mm; all panels use
# the same base font size, so that the text has the same size in all of them;
# the region PN is labelled PPN in the figure
figure_panels <- data.frame(panel = c('a', 'b', 'c', 'd', 'e', 'f', 'g'),
                            gene = c('RORA', 'ROBO1', 'LUC7L', 'FOXP2', 'FOXP2', 
                                     'GLUBK89', 'ENSGALG00010026609'),
                            side = 'L',
                            region = c('IMM', 'IMM', 'IMM', 'IMM', 'PN', 'IMM', 'IMM'),
                            width = c(half_width, half_width, 58, 58, 58, half_width, half_width),
                            height = c(78, 78, 58, 58, 58, 78, 78))
figure_pointsize <- 4.2 # largest text of StatPlot (cex = 2) at about 7 pt
figure_lwd_scale <- 0.4 # thinner lines than StatPlot's default (lwd = 2)

# full dataset
full_data_cols <- c('Name', 'TrUntr', 'Region', 'Side',
                    'Preference Score', 'Standardized Relative Value')
full_data <- data.frame(matrix(NA, 0, length(full_data_cols)))
colnames(full_data) <- full_data_cols

# inputs of the plots, for the panels of Figure 4
plot_inputs <- list()

# looping over folders
for(f in 1:length(folders)){

  # load data
  suppressMessages({
    behavioural <- as.data.frame(read_excel(file.path(data_folder, folders[f], 'behavioural.xls'), sheet = 1))
    predictors <- as.data.frame(read_excel(file.path(data_folder, folders[f], 'predictors.xls'), sheet = 1))
  })

  # harmonizing values
  colnames(predictors) <- c('Value', 'Sample')
  suppressWarnings({
    behavioural$Pref <- as.numeric(behavioural$Pref)
    predictors$Value <- as.numeric(predictors$Value)
  })

  # checking the data: leading or trailing spaces in the labels (e.g. 'Trained ')
  # would silently exclude samples from the groups
  behavioural[] <- lapply(behavioural, function(x) if(is.character(x)) trimws(x) else x)
  if(!all(behavioural$TrUntr %in% c('Trained', 'Untrained'))){
    stop('Unexpected values in TrUntr: ', folders[f])
  }
  if(!all(behavioural$Side %in% sides) || !all(behavioural$Region %in% regions)){
    stop('Unexpected values in Side or Region: ', folders[f])
  }

  # the two files must list the same samples in the same order
  if(!identical(behavioural$Sample, predictors$Sample)){
    stop('Not same samples: ', folders[f])
  }
  if(any(duplicated(behavioural$Sample))){
    message('Duplicated sample identifiers in ', folders[f], ': ',
            paste(unique(behavioural$Sample[duplicated(behavioural$Sample)]), collapse = ', '))
  }

  # ordering by sample
  sample_order <- order(behavioural$Sample)
  behavioural <- behavioural[sample_order, ]
  predictors <- predictors[sample_order, ]

  # x limits
  xLow <- floor(min(as.numeric(behavioural$Pref), na.rm = TRUE))
  xHigh <- 100 #ceiling(max(as.numeric(behavioural$Pref), na.rm = TRUE))

  # looping over brain parts
  for (s in sides){
    for(r in regions){

      # creating the results folder
      resFolder <- file.path(root_res_folder, folders[f],
                             paste0(ifelse(s == 'L', 'Left_', 'Right_'), r))
      dir.create(resFolder, showWarnings = FALSE, recursive = TRUE)

      # selecting the data
      idx <- which(behavioural$Side == s & behavioural$Region == r)
      behavioural_tmp <- behavioural[idx, ]
      predictors_tmp <- predictors[idx, ]

      # y limits
      yLow <- max(0, min(predictors_tmp$Value, na.rm = TRUE) - 0.1)
      yHigh <- max(predictors_tmp$Value, na.rm = TRUE) + 0.1

      # creating the dataset for the analysis
      dataset <- cbind(behavioural_tmp, value = predictors_tmp$Value)
      dataset <- dataset[!is.na(dataset$value), ]
      colnames(dataset)[ncol(dataset)] <- folders[f]

      # add to the full dataset
      tmp <- dataset[, c('TrUntr', 'Region', 'Side',
                         'Pref', folders[f])]
      tmp <- cbind(folders[f], tmp)
      colnames(tmp) <- full_data_cols
      full_data <- rbind(full_data, tmp)

      # inputs of the plot
      plot_inputs[[paste(folders[f], s, r)]] <- 
        list(dataset = dataset, 
             mainheading = paste(ifelse(s == 'L', 'Left', 'Right'), r), 
             ylabel = folders[f], responseno = dim(dataset)[2], 
             ylow = yLow, yhigh = yHigh, xlow = xLow, xhigh = xHigh)
      
      # analysis
      res <- tryCatch({
        png(filename = file.path(resFolder, 'figure.png'),
            width = 6000, height = 6000, res = 600)
        res <- StatPlot(dataset,
                        paste(ifelse(s == 'L', 'Left', 'Right'), r),
                        folders[f], dim(dataset)[2],
                        yLow, yHigh,
                        xLow, xHigh)
        dev.off()
        res
      },
      error = function(e){
        print(paste('Error in:', resFolder))
        dev.off()
        file.remove(file.path(resFolder, 'figure.png'))
        res <- NULL
      }
      )
      if(is.null(res)){
        next()
      }

      # writing results
      write.csv(res,
                file = file.path(resFolder, 'stats.csv'),
                row.names = FALSE)

    }

  }

}

# writing the full data
write.csv(full_data, row.names = FALSE,
          file = file.path(root_res_folder, 'full_data.csv'))

#### panels of Figure 4 ####

for(i in 1:nrow(figure_panels)){
  x <- figure_panels[i, ]
  inputs <- plot_inputs[[paste(x$gene, x$side, x$region)]]
  panel_folder <- paste0('Panel_', x$panel)
  dir.create(panel_folder, showWarnings = FALSE, recursive = TRUE)
  save_panel(file.path(panel_folder, paste0('panel_', x$panel)), 
             function() StatPlot(inputs$dataset, sub('PN$', 'PPN', inputs$mainheading), inputs$ylabel, 
                                 inputs$responseno, inputs$ylow, inputs$yhigh, 
                                 inputs$xlow, inputs$xhigh, 
                                 lwd_scale = figure_lwd_scale), 
             x$width, x$height, pointsize = figure_pointsize)
}
