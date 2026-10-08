#### Script for Figure 5f: analysis of the in situ hybridization data ####

# Figure 5f is the dotplot for IMM (Panel_f/IMM/dotplot_group.png); the same
# analysis is reported for NeoS (Panel_f/NeoS/)

##### set up #####

# memory and libraries
rm(list = ls())
library(tidyverse)
library(nlme)
library(emmeans)
source('../ancillary/figure_settings.R')
areas <- c('IMM', 'NeoS')

# results folder
res_folder <- 'Panel_f'
panel_width <- 58 # mm, three panels per row
panel_height <- 45 # mm

# reading the data (one value per cell)
dataset <- readRDS('../data/figure_5/totals.rds')

# summarize across cells
dataset <- dataset %>% group_by(group, side, area, batch, chick) %>%
                summarise(value = mean(value))

#### IMM ####
i <- 1

# settings
current_area <- areas[i]
current_res_folder <- file.path(res_folder, current_area)
current_dataset <- dataset %>% filter(area == current_area) %>%
  as.data.frame()
dir.create(current_res_folder, showWarnings = FALSE, recursive = TRUE)

# mixed model
# fixed factor: group, side
# random factor: batch / chick
m0 <- nlme::lme(value ~ group * side, random = ~ 1|batch/chick,
               data = current_dataset)

# side and interactions are not significant, let's remove them
m1 <- nlme::lme(value ~ group, random = ~ 1|batch/chick,
               data = current_dataset)

# recording results
sink(file.path(current_res_folder, 'results.txt'))

# analysis of variance of full model
cat('Full model')
print(anova(m0))
cat('\n\n')

# analysis of variance of final model
cat('Final model')
print(anova(m1))
cat('\n\n')

# difference between pairs
tmp <- emmeans(m1, 'group')
print(tmp)
cat('\n\n')
print(pairs(tmp))
cat('\n\n')

# stop sinking!
sink()

# dotplot by group
to_plot <- current_dataset %>% filter(area == current_area) %>%
  group_by(chick) %>% summarise(group = unique(group), value = mean(value))

p <- ggplot(data = to_plot,
            mapping = aes(x = group, y = value, color = group)) +
  geom_point(size = 1.5) +
  scale_x_discrete(name = 'Group', labels = c('Good learners', 'Poor learners', 'Untrained')) +
  scale_y_continuous(name = 'Normalized GLUBK89 Signal') +
  theme_bw(base_size = font_text) + theme_figure() +
  theme(legend.position = 'none') # the colors repeat the groups on the x axis
save_panel(file.path(current_res_folder, 'dotplot_group'), function() plot(p),
           panel_width, panel_height)

# writing the data of the plot
write.csv(to_plot, row.names = FALSE, file.path(current_res_folder, 'dotplot_values.csv'))

#### NeoS ####
i <- 2

# settings
current_area <- areas[i]
current_res_folder <- file.path(res_folder, current_area)
current_dataset <- dataset %>% filter(area == current_area) %>%
  as.data.frame()
dir.create(current_res_folder, showWarnings = FALSE, recursive = TRUE)

# mixed model
# fixed factor: group, side
# random factor: batch / chick
m <- nlme::lme(value ~ group * side, random = ~ 1|batch/chick,
               data = current_dataset)

# recording results
sink(file.path(current_res_folder, 'results.txt'))

# analysis of variance (no term is significant)
print(anova(m))
cat('\n\n')

# difference between pairs
tmp <- emmeans(m, 'group')
print(tmp)
cat('\n\n')
print(pairs(tmp))
cat('\n\n')

# stop sinking!
sink()

# dotplot by group
to_plot <- current_dataset %>% filter(area == current_area) %>%
  group_by(chick) %>% summarise(group = unique(group), value = mean(value))

p <- ggplot(data = to_plot,
            mapping = aes(x = group, y = value, color = group)) +
  geom_point(size = 1.5) +
  scale_x_discrete(name = 'Group', labels = c('Good learners', 'Poor learners', 'Untrained')) +
  scale_y_continuous(name = 'Normalized GLUBK89 Signal') +
  theme_bw(base_size = font_text) + theme_figure() +
  theme(legend.position = 'none') # the colors repeat the groups on the x axis
save_panel(file.path(current_res_folder, 'dotplot_group'), function() plot(p),
           panel_width, panel_height)

# writing the data of the plot
write.csv(to_plot, row.names = FALSE, file.path(current_res_folder, 'dotplot_values.csv'))
