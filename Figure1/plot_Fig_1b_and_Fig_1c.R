setwd("work_wd")
library(data.table)
library(dplyr)
library(ggplot2)
library(cowplot)

# read in ggplot theemes (made by Antonio Lentini) and plot parameters.
source("r_sources/AL_ggplot_themes.R")
source("r_sources/plot_parameters.R")

# read in data #
suppl_table_2 <- fread(file = "supplementary_tables/Supplementary_Table2.tsv")
suppl_table_3 <- fread(file = "supplementary_tables/Supplementary_Table3.tsv")

##### Figure 1b
# Calculate means
suppl_table_2_metrics <- suppl_table_2 %>% dplyr::group_by(individual, tissue_id) %>% rstatix::get_summary_stats(allelic_expression, type = "common")
suppl_table_2_metrics$tissue_id <- factor(suppl_table_2_metrics$tissue_id)
suppl_table_2_metrics <- suppl_table_2_metrics[order(suppl_table_2_metrics$mean, decreasing = F),]

###### Figure 1c ########
# calculate means
suppl_table_3_metrics <- suppl_table_3 %>% dplyr::group_by(individual, tissue_id) %>% rstatix::get_summary_stats(allelic_expression, type = "common")


ggsave(filename = paste0(plot_dir, "Fig1B_1C.pdf"), width = 100, height = 60, unit = "mm",
       plot_grid(ncol=2, rel_widths = c(1,0.35),
                 plot_grid(
                   ggplot(suppl_table_2_metrics, 
                          aes(x=reorder(individual, mean), y=mean)) + 
                     geom_ribbon(aes(ymin = mean-se, ymax = mean+se, group = 1), alpha = 0.15)+
                     geom_point(aes(col = tissue_id, fill = tissue_id), size = 0.75) + 
                     geom_hline(yintercept = 0.4, lty = 2) + 
                     theme_AL_box_rotX() + 
                     theme(legend.position = "none", 
                           legend.title = element_blank(), 
                           axis.text.x = element_text(size=6), 
                           axis.text.y = element_text(size=6), 
                           axis.title.x = element_text(size=6), 
                           axis.title.y = element_text(size=6))+
                     labs(x="participant", y="allelic expression (mean)")+
                     scale_y_continuous(breaks=c(0,0.1,0.2,0.3,0.4,0.5), limits = c(0,0.5))+
                     scale_color_manual(values=color_valzz, limits = force)+
                     geom_text(data=suppl_table_2_metrics[suppl_table_2_metrics$individual %in% suppl_table_3_metrics$individual,], 
                               aes(x=individual, y=0.3, label = individual), angle = 90, size = 2)),
                 
                 plot_grid(ggplot(suppl_table_3_metrics, aes(x=reorder(individual, mean), y=mean)) + 
                             stat_summary(geom="crossbar",
                                          fun.min = function(z) {quantile(z,0.25)},
                                          fun.max = function(z) {quantile(z,0.75)},
                                          fun = median,
                                          width = .75) + 
                             geom_hline(yintercept = 0.4, lty = 2) + 
                             theme_AL_box_rotX() + 
                             theme(axis.text.x = element_text(size=6), 
                                   axis.text.y = element_text(size=6), 
                                   axis.title.x = element_text(size=6), 
                                   axis.title.y = element_text(size=6))+
                             scale_y_continuous(breaks=c(0,0.1,0.2,0.3,0.4,0.5), limits = c(0,0.5))+
                             labs(x="participant", y="allelic expression (mean)"))
       )
)




