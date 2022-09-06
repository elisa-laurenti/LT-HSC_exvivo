suppressMessages(library(ggplot2))

meta_plot = read.csv('meta_PD_QA81_22_N_mPB_CycleRegress__UMAP_forPlotting.csv',row.names=1)

meta_pseudo = read.csv('meta_PD_QA81_22_N_mPB__CycleRegress_pseudotime_monocle3.csv',row.names=1)

identical(rownames(meta_plot),rownames(meta_pseudo))

meta_pseudo$timepoints = meta_plot$both_time
table(meta_pseudo$timepoints)

color_dir=c('LT_0h'='blue','LT_6h'='cornflowerblue',
             'LT_24h_UNTR'= 'mediumseagreen', 
             'LT_72h_UNTR'= 'darkred',
             'LT_24h_PD'='blueviolet',
             'LT_72h_PD'='hotpink',
             'Day0_LT-HSC_NT'='darkorange',
             'Day3_LT-HSC_GFP+PD-'='lawngreen',
             'Day3_LT-HSC_GFP+PD+'='black')
options(repr.plot.width=11, repr.plot.height=8)

ggplot(meta_pseudo, aes(x=pseudotime_2D_rank, color=timepoints)) +
  geom_density()+ggtitle('PD_QA81_22_N_mPB__CycleRegress pseudotime 2D ranking density')+ 
  scale_color_manual(values = color_dir)


