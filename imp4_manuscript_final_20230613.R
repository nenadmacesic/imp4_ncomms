#Load relevant packages
library(tidyverse)
library(reshape2)
library(lubridate)
library(ape)
library(phangorn)
library(ggtree)
library(ggtreeExtra)
library(ggnewscale)
library(RColorBrewer)
library(treeio)
library(phytools)
library(scales)
library(ggridges)
library(zoo)
library(aplot)
library(svglite)
library(patchwork)
library(igraph)
library(ggraph)
library(ggnewscale)
library(ggforce)

setwd('/your_directory')

#*Load metadata####
metadata_final <- read_csv('metadata_final_20230613.csv')

#*Figure 1A####
pie_chart_df <- group_by(metadata_final, species_st_final) %>%
                summarise(total_iso_final = n())

species_st_pie_chart <- ggplot(pie_chart_df, aes(x="", y=total_iso_final, fill=species_st_final))+
                       geom_bar(width = 1, stat = "identity") +
                       coord_polar("y", start=0) +
                       theme_minimal() +
                       scale_fill_brewer(palette= "Paired") +
                       theme(axis.title = element_blank(),
                              axis.text = element_blank(),
                              title = element_text(size = 14),
                              legend.title = element_text(size = 12, face = 'bold'),
                              #legend.position = 'bottom',
                              legend.text =  element_text(size = 12)) +
                        #ggtitle('Bacterial host strain proportion') +
                        labs(fill = 'Bacterial host strain')

#*Figure 1B####
metadata_year_species_summary <-  group_by(metadata_final, species_st_clean, rep_type_clean) %>%
                                  mutate(total_isos = n()) %>%
                                  distinct(species_st_clean, rep_type_clean, .keep_all = TRUE) %>%
                                  select(species_clean, species_st_clean, rep_type_clean, total_isos)

metadata_year_species_summary_period <- group_by(metadata_final, species_st_clean, rep_type_clean, period) %>%
                                        mutate(total_isos = n()) %>%
                                        distinct(species_st_clean, rep_type_clean, period, .keep_all = TRUE) %>%
                                        select(species_st, rep_type_clean, species_clean, species_st_clean, total_isos, period)

# Define colour palette and number of species for plotting
nb.cols_source <- length(unique(metadata_final$rep_type_clean))
colors_source <- colorRampPalette(brewer.pal(nb.cols_source, "Set1"))(nb.cols_source)
fill_source <- colorRampPalette(brewer.pal(nb.cols_source, "Set1"))(nb.cols_source)

nb.cols_source <- length(unique(metadata_final$rep_type_clean)) ## changed to var here
fill_source <- colorRampPalette(brewer.pal(nb.cols_source, "Set1"))(nb.cols_source) %>%
                as.data.frame() %>%
                rename("fill_colour" = ".") %>%
                add_column("Genetic setting" = sort(unique(metadata_final$rep_type_clean)))

colors_source <- colorRampPalette(brewer.pal(nb.cols_source, "Set1"))(nb.cols_source) %>%
                  as.data.frame() %>%
                  rename("fill_colour" = ".") %>%
                  add_column("Genetic setting" = sort(unique(metadata_final$rep_type_clean)))

species_st_year_barplot <- ggplot(data=metadata_year_species_summary, aes(x=species_st_clean, y=total_isos, fill = rep_type_clean)) +
                            geom_bar(stat="identity") +
                            scale_fill_brewer(palette= "Set1") +
                            facet_grid(vars(species_clean), scales = 'free', space = 'free') +
                            theme(strip.text.y = element_blank(), 
                                  axis.text = element_text(size = 14),
                                  legend.text = element_text(size = 14)) +
                            coord_flip()

species_st_year_barplot_period <- ggplot(data=metadata_year_species_summary_period, aes(x=species_st_clean, y=total_isos, fill = rep_type_clean)) +
                                  geom_bar(stat="identity") +
                                  scale_fill_manual("Genetic setting", values = setNames(as.character(fill_source$fill_colour),fill_source$`Genetic setting`), drop = TRUE) +
                                  scale_y_continuous(breaks = pretty_breaks()) +
                                  facet_grid(vars(species_clean), vars(period), scales = 'free', space = 'free') +
                                  theme_bw() +
                                  theme(strip.text.y = element_blank(), 
                                        strip.text.x = element_text(size = 16, face = 'bold'),
                                        axis.text = element_text(size = 16),
                                        axis.title = element_text(size = 16, face = "bold"),
                                        legend.text = element_text(size = 16),
                                        legend.title = element_text(size = 16, face = "bold")) +
                                  ylab("Total no. genomes") +
                                  xlab("Bacterial host strain") +
                                  coord_flip()

ggsave('species_st_year_barplot_period_20221205.svg', plot=species_st_year_barplot_period,
       width=6000, height=3750, units = "px")

#*Figure 2A####
serratia_alfred_tree <- read.tree('imp4_serratia_phylo_20210928.newick')
serratia_tip_labels_to_drop <- (serratia_alfred_tree$tip.label)[str_detect(serratia_alfred_tree$tip.label, "GCF")]
serratia_alfred_tree_2 <- drop.tip(serratia_alfred_tree, serratia_tip_labels_to_drop)

serratia_alfred_ggtree <- ggtree(serratia_alfred_tree_2, #layout = 'circular', branch.length = 'none',
                                 color = "black", size = 0.5, right = T) +
                          geom_tiplab(color = NA, align=TRUE, offset=.0001) + # can remove this if   tip dots not required - useful for heatmaps though
                          geom_cladelab(node = 100, label = "Lineage 2", align = TRUE, angle = 90, hjust = 'center', offset.text = 0.0005) +         
                          geom_cladelab(node = 54, label = "Lineage 1", align = TRUE, angle = 90, hjust = 'center', offset.text = 0.0005)         

serratia_metadata <- as.data.frame(get_taxa_name(serratia_alfred_ggtree))
colnames(serratia_metadata) <- 'genome'
serratia_metadata <- mutate(serratia_metadata, isolate = str_extract(genome, "CPO[0-9]{3}"),
                            isolate = ifelse(genome=='Reference', "CPO051", isolate)) %>%
                     left_join(., metadata_final, by = 'isolate') %>%
                     mutate(environmental = case_when(environmental == '1' ~ 'Environmental',
                                                     environmental == '0' ~ 'Clinical')) # convert to factored name for ease of plotting etc          

serratia_year_long <- select(serratia_metadata, genome, year) %>%
                      mutate("variable" = "Year")

# Produce heatmap for year
serratia_year_hmap <- ggplot(serratia_year_long, aes(x=as.factor(variable), y=genome)) + 
                      geom_tile(fill = "white", color = "grey") + # width = 0.1
                      geom_text(aes(label = year), size = 3) +
                      theme_tree2(axis.text.x = element_text(angle = 45, hjust=1, size = 9)) +
                      coord_fixed(ratio = 0.2) # this sets the total width of the heatmap, nice.

# Prepare mobtype_cluster df
serratia_rep_long <- select(serratia_metadata, genome, rep_type_clean) %>%
                     mutate("variable" = "Genetic setting")

# Produce heatmap for mobtype_cluster
serratia_rep_hmap <- ggplot(serratia_rep_long, aes(x=as.factor(variable), y=genome)) + 
                      geom_tile(aes(fill=as.factor(rep_type_clean)), color = "grey") +
                      scale_fill_manual("Genetic setting", values = setNames(as.character(fill_source$fill_colour),fill_source$`Genetic setting`), drop = TRUE) +
                      theme_tree2(axis.text.x = element_text(angle = 45, hjust=1, size = 9))

serratia_enviro_long <- select(serratia_metadata, genome, environmental) %>%
                        mutate('variable' = 'Source')

serratia_enviro_hmap <- ggplot(serratia_enviro_long, aes(x=as.factor(variable), y=genome)) + 
                        geom_tile(aes(fill=as.factor(environmental)), color = "grey") +
                        scale_fill_brewer(name = "Source", palette = "Greys") +
                        theme_tree2(axis.text.x = element_text(angle = 45, hjust=1, size = 9))

serratia_combined <- serratia_year_hmap %>%
                      insert_right(serratia_rep_hmap, width = 0.15) %>%
                      insert_right(serratia_enviro_hmap, width = 0.15) %>%
                      insert_left(serratia_alfred_ggtree)

ggsave(plot = serratia_combined, file = 'alfred_serratia_tree_big_20220524.svg', width = 2000, height = 4000, unit = 'px')

#*Figure 2A inset####
chromo_tree_fn <- function(x, y, z){
  #Input: x = phylogenetic tree
  #Input: y = size of points
  #Input: z = figure title
  alfred_tree <- read.tree(x)
  
  alfred_tree_genomes <- as.data.frame(alfred_tree$tip.label)
  colnames(alfred_tree_genomes) <- "genome"
  alfred_tree_genomes <- mutate(alfred_tree_genomes,
                                origin = ifelse(str_detect(genome, "^CPO"), "Alfred", "Non-Alfred"),
                                origin = ifelse(genome=="Reference", "Alfred", origin))
  
  alfred_ggtree <- ggtree(alfred_tree,
                          color = "black", size = 0.5, right = T) +
    theme(legend.position = 'none')
  
  alfred_ggtree_2 <- alfred_ggtree %<+% alfred_tree_genomes +
    geom_tippoint(aes(colour = origin), size = y) +
    geom_treescale(y = -2) +
    theme(title = element_text(size = 14, face = 'bold', hjust = 0.5)) +
    ggtitle(z)
}

serratia_tree <- chromo_tree_fn('imp4_serratia_phylo_20210928.newick', 1, expression(italic("")))

#*Figure 2B####
mobtyper_results_imp4_contigs_mashtree_plot <- group_by(metadata_final, species_st) %>%
                                                mutate(total_species_st = n(),
                                                       species_st_clean = paste(species_final, " ", "ST", st, sep = ""),
                                                       species_st_clean = ifelse(total_species_st>2, species_st_clean, paste(species_clean, " - Other STs", sep = "")),
                                                       species_st_clean = ifelse(str_detect(species_st_clean, "Acinetobacter"), "Acinetobacter spp.", species_st_clean),
                                                       species_st_clean = ifelse(str_detect(species_st_clean, "Citrobacter"), "Citrobacter spp.", species_st_clean),
                                                       species_st_clean = ifelse(str_detect(species_st_clean, "Serratia marcescens ST1"), "Serratia marcescens lin. 1", species_st_clean),
                                                       species_st_clean = ifelse(str_detect(species_st_clean, "Serratia marcescens ST2"), "Serratia marcescens lin. 2", species_st_clean)) %>%
                                                ungroup()

mobtyper_results_imp4_contigs_mashtree_plot_2 <- distinct(mobtyper_results_imp4_contigs_mashtree_plot,
                                                          species_st, .keep_all = TRUE) %>%
                                                 select(species_st, species_st_clean)

### Input trees
incC_tree_ben <- read.tree('imp4_incC_no_fusion_mashtree_20221222.newick')
incC_tree_ben$tip.label <- str_extract(incC_tree_ben$tip.label, "CPO[0-9]{3}")

incHI2_tree_ben <- read.tree('imp4_incHI2_no_fusion_mashtree_20221222.newick')
incHI2_tree_ben$tip.label <- str_extract(incHI2_tree_ben$tip.label, "CPO[0-9]{3}")

incC_tree_r <- phytools::midpoint.root(incC_tree_ben)
incHI2_tree_r <- phytools::midpoint.root(incHI2_tree_ben)

tip_list_ordered_incHI2 <- as.data.frame(incHI2_tree_r$tip.label) %>% # produce tip list with ape and output as df
  rename("tip_names" = "incHI2_tree_r$tip.label")

tip_list_ordered_incC <- as.data.frame(incC_tree_r$tip.label) %>% # produce tip list with ape and output as df
  rename("tip_names" = "incC_tree_r$tip.label")

tip_list_ordered_incHI2_incC <- bind_rows(tip_list_ordered_incHI2, tip_list_ordered_incC)

metadata_incHI2_incC <-  filter(metadata_final, isolate %in% tip_list_ordered_incHI2_incC$tip_names) %>% # Filter data from df that aren't in the phylogeny tips
                         select(isolate, species_st, species_st_clean, rep_type_clean, flanker_cluster, nt_seq, environmental, year) %>%
                         mutate(species_st_clean = ifelse(species_st_clean=='Klebsiella pneumoniae ST340', "Klebsiella pneumoniae complex - Other STs",
                                                          species_st_clean),
                                species_st_clean = ifelse(species_st_clean=='Escherichia coli ST176', "Escherichia coli - Other STs",
                                                          species_st_clean),
                                environmental = case_when(environmental == '1' ~ 'Environmental',
                                                          environmental == '0' ~ 'Clinical')) # convert to factored name for ease of plotting etc

metadata_incHI2_incC_flanker <- select(metadata_incHI2_incC, isolate, flanker_cluster) %>%
                                mutate(presence = 1,
                                       variable = 'Flanking cluster') %>%
                                dcast(isolate ~ flanker_cluster, value.var = 'presence') %>%
                                melt(., id.var = 'isolate') %>%
                                mutate(value = replace_na(value, 0),
                                       value = ifelse(value==0, "Absent", 'Present')) %>%
                                rename(flanker_cluster = variable, presence = value)

nb.cols_species_st_incHI2_incC <- length(unique(metadata_incHI2_incC$species_st_clean)) # change 'species_final' to other variable to get new list.
fill_species_st_incHI2_incC <- colorRampPalette(brewer.pal(4, "Set1"))(nb.cols_species_st_incHI2_incC) %>%
                                as.data.frame() %>%
                                rename("fill_colour" = ".") %>%
                                add_column("Bacterial strain" = sort(unique(metadata_incHI2_incC$species_st_clean)))

nb.cols_nt_seq_incHI2_incC <- length(unique(metadata_incHI2_incC$nt_seq)) # change 'species_final' to other variable to get new list.
fill_nt_seq_incHI2_incC <- colorRampPalette(brewer.pal(8, "Dark2"))(nb.cols_nt_seq_incHI2_incC) %>%
                           as.data.frame() %>%
                           rename("fill_colour" = ".") %>%
                           add_column("Integron SNV group" = sort(unique(metadata_incHI2_incC$nt_seq)))

mashtree_function_nenad <- function(x, y){
  #Input: 'x' is Ben's mashtree edit
  #Input: 'y' is rep type that want to focus on
  
  
  # Construct tree with species as tips
  mdata_merge_incHI2 <- filter(metadata_incHI2_incC, str_detect(rep_type_clean, y))  
  
  incHI2_tree_r_mdata <- ggtree(x) %<+% mdata_merge_incHI2 + # attach metadata like this
    geom_tippoint(aes(fill = species_st_clean,
                      colour = species_st_clean), size = 5) +
    #shape = as.factor(environmental),
    #size=5)) +
    #scale_shape_manual("Isolate type", values = c(21, 22)) +
    scale_fill_manual("Bacterial strain", values = setNames(as.character(fill_species_st_incHI2_incC$fill_colour),
                                                            fill_species_st_incHI2_incC$`Bacterial strain`), drop = TRUE) +
    scale_color_manual("Bacterial strain", values = setNames(as.character(fill_species_st_incHI2_incC$fill_colour),
                                                             fill_species_st_incHI2_incC$`Bacterial strain`), drop = TRUE) +
    guides(fill = guide_legend(override.aes = list(shape = 21)),
           shape = guide_legend(override.aes = list(fill = "black"))
    ) +
    #geom_tiplab(aes(color = species_st), align=TRUE)
    geom_tiplab(color = NA, align=TRUE, offset=.0001) + # can remove this if   tip dots not required - useful for heatmaps though
    #geom_treescale(fontsize = 4, x = 0.001, y = 40) + # can remove this if   scale not required
    #geom_treescale(fontsize = 4) + # can remove this if   scale not required
    theme_tree()
  
  ### Prepare data for heatmap presentation for second plot
  
  # Prepare year df
  year_long <- mdata_merge_incHI2 %>%
    select(isolate, year) %>%
    mutate("variable" = "Year")
  
  # Produce heatmap for year
  year_hmap <- ggplot(year_long, aes(x=as.factor(variable), y=isolate)) +
    geom_tile(fill = "white", color = "grey") + # width = 0.1
    geom_text(aes(label = year), size = 3) +
    theme_tree2(axis.text.x = element_text(angle = 45, hjust=1, size = 9)) +
    coord_fixed(ratio = 0.2) # this sets the total width of the heatmap, nice.
  
  #Produce heatmap for environmental
  enviro_long <- select(mdata_merge_incHI2, isolate, environmental) %>%
    mutate('variable' = 'Source')
  
  enviro_hmap <- ggplot(enviro_long, aes(x=as.factor(variable), y=isolate)) + 
    geom_tile(aes(fill=as.factor(environmental)), color = "grey") +
    scale_fill_brewer(name = "Source", palette = "Greys") +
    theme_tree2(axis.text.x = element_text(angle = 45, hjust=1, size = 9))
  
  # Prepare mobtype_cluster df
  mob_long <- mdata_merge_incHI2 %>%
    select(isolate, rep_type_clean) %>%
    mutate("variable" = "Genetic setting")
  
  # Produce heatmap for mobtype_cluster
  mob_hmap <- ggplot(mob_long, aes(x=as.factor(variable), y=isolate)) +
    geom_tile(aes(fill=as.factor(rep_type_clean)), color = "grey") +
    scale_fill_manual("Genetic setting", values = setNames(as.character(fill_source$fill_colour),
                                                           fill_source$`Genetic setting`), drop = TRUE) +
    # scale_fill_brewer(name = "Mobtype cluster", palette = "RdPu", trans =   "reverse") +
    theme_tree2(axis.text.x = element_text(angle = 45, hjust=1, size = 9))

  flank_long <- filter(metadata_incHI2_incC_flanker, isolate %in% mdata_merge_incHI2$isolate)
  
  flank_hmap <- ggplot(flank_long, aes(x=as.factor(flanker_cluster), y=isolate)) +
    geom_tile(aes(fill=as.factor(presence)), colour = "grey") +
    scale_fill_manual(name = "Flanking cluster", values = c("#FFFFFF", "#000000")) + # for black and white fill
    scale_x_discrete("Flanking cluster") +
    theme_tree2(axis.text.x = element_text(angle = 45, hjust=1, size = 9)) +
    theme(axis.text.x = element_blank(),
          axis.ticks = element_blank(),
          axis.title.x = element_text(angle = 45, hjust=1, size = 9))
  
  # Prepare nt_seq_group df
  nt_seq_long <- mdata_merge_incHI2 %>%
    select(isolate, nt_seq) %>%
    mutate("variable" = "Integron SNV group")
  
  # Produce heatmap for nt_seq_group
  nt_seq_hmap <- ggplot(nt_seq_long, aes(x=as.factor(variable), y=isolate)) +
    geom_tile(aes(fill=as.factor(nt_seq)), color = "grey") +
    scale_fill_manual("Integron SNV group", values = setNames(as.character(fill_nt_seq_incHI2_incC$fill_colour),
                                                              fill_nt_seq_incHI2_incC$`Integron SNV group`), drop = TRUE) +
    theme_tree2(axis.text.x = element_text(angle = 45, hjust=1, size = 9))
  
  # Combine all with aplot
  combined <- year_hmap %>%
    insert_right(enviro_hmap, width = 0.15) %>%
    insert_right(mob_hmap, width = 0.15) %>% # only issue here is my
    insert_right(flank_hmap, width = 0.35) %>%
    insert_right(nt_seq_hmap, width = 0.15) %>%
    insert_left(incHI2_tree_r_mdata)
  combined
}

incC_mashtree_nenad <- mashtree_function_nenad(incC_tree_r, "IncC")

ggsave(file="incC_mashtree_nenad_20221223.svg", plot=incC_mashtree_nenad, width=2500, height=4000, units = 'px')

#*Figure 3A####
alfred_tree_fn <- function(x, y, z){
  #Input: 'x' is Newick tree to draw
  #Input: 'y' is tip points to specifically exclude
  #Input: 'z' is name of Reference isolate
  
  alfred_tree <- read.tree(x)
  tip_labels_to_drop <- c((alfred_tree$tip.label)[str_detect(alfred_tree$tip.label, "GCF")], y)
  
  alfred_tree_2 <- drop.tip(alfred_tree, tip_labels_to_drop)
  
  alfred_ggtree <- ggtree(alfred_tree_2, #layout = 'circular', branch.length = 'none',
                          color = "black", size = 0.5, right = T) +
    geom_tiplab(color = NA, align=TRUE, offset=.0001) #+# can remove this if   tip dots not required - useful for heatmaps though
  #geom_treescale(x = 0, y = -10, offset = 1)
  
  trial_metadata <- as.data.frame(get_taxa_name(alfred_ggtree))
  colnames(trial_metadata) <- 'genome'
  trial_metadata <- mutate(trial_metadata, isolate = str_extract(genome, "CPO[0-9]{3}"),
                           isolate = ifelse(genome=='Reference', z, isolate)) %>%
    left_join(., metadata_final, by = 'isolate') %>%
    mutate(environmental = case_when(environmental == '1' ~ 'Environmental',
                                     environmental == '0' ~ 'Clinical')) # convert to factored name for ease of plotting etc          
  
  year_long <- select(trial_metadata, genome, year) %>%
    mutate("variable" = "Year")
  
  # Produce heatmap for year
  year_hmap <- ggplot(year_long, aes(x=as.factor(variable), y=genome)) + 
    geom_tile(fill = "white", color = "grey") + # width = 0.1
    geom_text(aes(label = year), size = 3) +
    theme_tree2(axis.text.x = element_text(angle = 45, hjust=1, size = 9)) +
    coord_fixed(ratio = 0.2) # this sets the total width of the heatmap, nice.
  
  # Prepare mobtype_cluster df
  rep_long <- select(trial_metadata, genome, rep_type_clean) %>%
    mutate("variable" = "Genetic setting")
  
  # Produce heatmap for mobtype_cluster
  rep_hmap <- ggplot(rep_long, aes(x=as.factor(variable), y=genome)) + 
    geom_tile(aes(fill=as.factor(rep_type_clean)), color = "grey") +
    scale_fill_manual("Genetic setting", values = setNames(as.character(fill_source$fill_colour),fill_source$`Genetic setting`), drop = TRUE) +
    theme_tree2(axis.text.x = element_text(angle = 45, hjust=1, size = 9)) +
    theme(legend.position = 'none')
  
  enviro_long <- select(trial_metadata, genome, environmental) %>%
    mutate('variable' = 'Source')
  
  enviro_hmap <- ggplot(enviro_long, aes(x=as.factor(variable), y=genome)) + 
    geom_tile(aes(fill=as.factor(environmental)), color = "grey") +
    scale_fill_brewer(name = "Source", palette = "Greys") +
    theme_tree2(axis.text.x = element_text(angle = 45, hjust=1, size = 9)) +
    theme(legend.position = 'none')
  
  combined <- year_hmap %>%
    insert_right(rep_hmap, width = 0.15) %>%
    insert_right(enviro_hmap, width = 0.15) %>%
    insert_left(alfred_ggtree)
}

alfred_st114_tree_big <- alfred_tree_fn('entcpx114.clean.core.aln.treefile', "", "CPO009")
ggsave(plot = alfred_st114_tree_big, file = 'alfred_st114_tree_big_20220524.svg', width = 1000, height = 1500, units = 'px')

#*Figure 3A inset####
alfred_st114_tree <- chromo_tree_fn('entcpx114.clean.core.aln.treefile', 1, expression(paste(italic(""), "")))
ggsave(plot = alfred_st114_tree, file = 'alfred_st114_tree_20220524.svg', height = 1000, width = 500, unit = 'px')

#*Figure 3B####
flanker_cluster_tree_df <- read_csv('flanker_output_2_edit_20221219_clean.csv') %>%
                           mutate(nenad_cluster = as.factor(nenad_cluster),
                                  flanker_window = as.integer(flanker_window))

flanker_output <- read_csv('flanker_output_20230613.csv')

flanker_cluster_pie_chart_df <- select(flanker_output, isolate, flanker_window_cluster, flanker_window,
                                       flanker_cluster) %>%
                                filter(!is.na(flanker_cluster),
                                       flanker_window %in% c(5000)) %>%
                                left_join(., select(metadata_final, isolate, rep_type_clean), by = 'isolate') %>%
                                group_by(flanker_window_cluster, rep_type_clean) %>%
                                summarise(total_contigs = n()) %>%
                                group_by(flanker_window_cluster) %>%
                                mutate(total_cluster = sum(total_contigs)) %>%
                                filter(total_cluster>4) %>%
                                left_join(., select(flanker_cluster_tree_df, flanker_window_cluster, nenad_cluster),
                                          by = 'flanker_window_cluster')

cp <- coord_polar(theta = "y")
cp$is_free <- function() TRUE

flanker_cluster_pie_chart <- ggplot(flanker_cluster_pie_chart_df, aes(x="", y = total_contigs, fill=rep_type_clean))+
                            geom_bar(width = 1, stat = "identity") +
                            cp +
                            facet_wrap(~flanker_window_cluster, ncol=4, nrow=4, scales = 'free') +
                            theme_minimal() +
                            scale_fill_manual("Genetic setting", values = setNames(as.character(fill_source$fill_colour),
                                                                                   fill_source$`Genetic setting`), drop = TRUE) +                    
                            theme(aspect.ratio = 1,
                                  axis.title = element_blank(),
                                  axis.text = element_blank(),
                                  title = element_text(size = 14),
                                  legend.title = element_text(size = 12, face = 'bold'),
                                  legend.text =  element_text(size = 12)) +
                            labs(fill = 'Genetic setting')

ggsave(plot = flanker_cluster_pie_chart, file = 'flanker_cluster_pie_chart_20221220.svg', width = 5000, height = 4000, unit = 'px')

#*Figure 4A####
incHI2_mashtree_nenad <- mashtree_function_nenad(incHI2_tree_r, "IncHI2A")

ggsave(file="incHI2_mashtree_nenad_20221223.svg", plot=incHI2_mashtree_nenad, width=2500, height=4000, units = 'px')

#*Figure 5A####
mobtyper_results_imp4_contigs_4_pts <- filter(metadata_final, !is.na(ur_dummy) & environmental==0) %>%
                                        group_by(ur_dummy) %>%
                                        mutate(acc_coll_date_min = min(acc_coll_date),
                                               time_from_first_iso = acc_coll_date-acc_coll_date_min)

mobtyper_results_imp4_contigs_4_pts_2 <- ungroup(mobtyper_results_imp4_contigs_4_pts) %>%
                                          select(ur_dummy, isolate, species_st, time_from_first_iso, species_st_clean, rep_type_clean) %>%
                                          group_by(ur_dummy) %>%
                                          mutate(time_from_first_iso = as.integer(time_from_first_iso),
                                                 total_rep_type = n_distinct(rep_type_clean),
                                                 total_species_st = n_distinct(species_st),
                                                 patient_class = ifelse(total_rep_type>1,
                                                                        "Multiple col.\n events", NA),
                                                 patient_class = ifelse(total_rep_type==1&total_species_st>1,
                                                                        "Inter-strain plasmid transfer", patient_class),
                                                 patient_class = ifelse(is.na(patient_class), "Persistent colonization", patient_class)) %>%
                                          group_by(ur_dummy, time_from_first_iso) %>%
                                          mutate(overlap_iso = n(),
                                                 overlap_iso_no = row_number(),
                                                 time_from_first_iso_2 = ifelse(overlap_iso>1&overlap_iso_no == 2, 2, time_from_first_iso),
                                                 time_from_first_iso_2 = ifelse(overlap_iso>1&overlap_iso_no == 3, 4, time_from_first_iso_2),
                                                 diff_iso = ifelse(time_from_first_iso_2==time_from_first_iso, 0, 1))


shape_rep_type_clean <- distinct(mobtyper_results_imp4_contigs_4_pts_2, rep_type_clean) %>%
                        mutate(shape = rep_type_clean,
                               shape = case_when(shape=='Chromosome'~'19',
                                                 shape=='IncC'~'15',
                                                 shape=='IncFIA/IncFIB/IncP'~'4',
                                                 shape=='IncHI2A type 1'~'17',
                                                 shape=='IncHI2A type 2'~'18',
                                                 shape=='IncL/M'~'16'),
                               shape = as.integer(shape))

nb.cols_species_st_within_pt <- length(unique(mobtyper_results_imp4_contigs_4_pts_2$species_st_clean)) # change 'species_final' to other variable to get new list.
fill_species_st_within_pt <- colorRampPalette(brewer.pal(8, "Set1"))(nb.cols_species_st_within_pt) %>%
                              as.data.frame() %>%
                              rename("fill_colour" = ".") %>%
                              add_column("Species/ST" = sort(unique(mobtyper_results_imp4_contigs_4_pts_2$species_st_clean)))



within_pt_plot <- ggplot(mobtyper_results_imp4_contigs_4_pts_2, aes(x=time_from_first_iso_2, y=ur_dummy, shape = rep_type_clean,
                                                                    color = species_st_clean)) + 
                  geom_point(size = 3) +
                  scale_color_manual("Bacterial strain", values = setNames(as.character(fill_species_st_within_pt$fill_colour),
                                                                           fill_species_st_within_pt$`Species/ST`), drop = TRUE) +
                  scale_shape_manual("Genetic setting", values = setNames(shape_rep_type_clean$shape,
                                                                          shape_rep_type_clean$rep_type_clean), drop = TRUE) +
                  facet_grid(rows = vars(patient_class), scales = 'free', space = 'free') +
                  theme(axis.title=element_text(face='bold'),
                        legend.title = element_text(face = 'bold'),
                        strip.text = element_text(face = 'bold')) +
                  xlab("Time from first blaIMP-4 isolate") +
                  ylab("Patient ID") +
                  labs(color = 'Bacterial strain',
                       shape = "Genetic setting") 

ggsave(plot = within_pt_plot, file = 'within_pt_plot_20221223.svg', height = 3000, width = 3750, unit = 'px')

#*Figure 5B####
non_dup_fn_pt_mvmt <- function(x){
  x [1:2] <- t(apply(x[1:2], 1, sort))
  snp_matrix_unique <- distinct(x, ur_dummy.x, ur_dummy.y, .keep_all = TRUE)  
  return(snp_matrix_unique)
}

seq_plan_metadata_imp4_2_non_env <- filter(metadata_final, environmental==0)

seq_plan_metadata_imp4_2_non_env_2 <- filter(seq_plan_metadata_imp4_2_non_env,
                                             !is.na(ur_dummy))

pt_first_iso <- group_by(seq_plan_metadata_imp4_2_non_env, ur_dummy) %>%
                mutate(first_iso = min(acc_coll_date)) %>%
                filter(acc_coll_date==first_iso) %>%
                select(ur_dummy, first_iso_date = acc_coll_date)


#**SNP distances
snp_distance_fn <- function(x){
  snp_distances <- read_csv(x)
  snp_distances_melt <- melt(snp_distances)
  colnames(snp_distances_melt) <- c('isolate_1', 'isolate_2', 'value')
  
  snp_distances_melt_2 <- filter(snp_distances_melt, !(isolate_1==isolate_2))
  snp_distances_melt_2[1:2] <- t(apply(snp_distances_melt_2[1:2], 1, sort))
  snp_distances_melt_3 <- distinct(snp_distances_melt_2, isolate_1, isolate_2, .keep_all = TRUE)
  snp_distances_melt_4 <- mutate(snp_distances_melt_3,
                                 isolate_1_class = ifelse(str_detect(isolate_1, "^GCF"), 2, 1),
                                 isolate_2_class = ifelse(str_detect(isolate_2, "^GCF"), 2, 1),
                                 isolate_rel = isolate_1_class+isolate_2_class,
                                 isolate_rel = as.character(isolate_rel),
                                 isolate_rel = ifelse(isolate_rel=='2', "Alfred", isolate_rel),
                                 isolate_rel = ifelse(isolate_rel=='3', "Mixed", isolate_rel),
                                 isolate_rel = ifelse(isolate_rel=='4', "Non-Alfred", isolate_rel),
                                 isolate_rel = as.factor(isolate_rel))
  
  snp_distances_summary <- group_by(snp_distances_melt_4, isolate_rel) %>%
    summarize(med_distance = as.vector(summary(value)[3]),
              mean_distance = as.vector(summary(value)[4]),
              first_quart_distance = as.vector(summary(value)[2]),
              third_quart_distance = as.vector(summary(value)[5])) %>%
    mutate(origin = x)
  snp_matrix_joint <- list(snp_distances_melt_4, snp_distances_summary)
  
}

snp_dist_clean_fn <- function(x, y){
  #Input: x = CSV of SNP distance matrix
  #       y = Isolate code of reference
  snp_distances <- (snp_distance_fn(x))[[1]]
  snp_distances <- rename(snp_distances, isolate.x = isolate_1, isolate.y = isolate_2) %>%
    filter(isolate_rel == "Alfred") %>%
    mutate(isolate.x = str_extract(isolate.x, "CPO[0-9]{3}"),
           isolate.y = str_extract(isolate.y, "CPO[0-9]{3}"),
           isolate.y = ifelse(is.na(isolate.y), y, isolate.y)) %>%
    select(isolate.x, isolate.y, snp_dist = value)
}

snp_distances_serratia_lin_1_clean <- snp_dist_clean_fn('snp_matrix_serratia_lin_1.csv', 'CPO101')
snp_distances_serratia_lin_2_clean <- snp_dist_clean_fn('snp_matrix_serratia_lin_2.csv', 'CPO037')
snp_distances_pseudomonas_clean <- snp_dist_clean_fn('snp_matrix_pseudomonas_st111.csv', 'CPO100')
snp_distances_entcpx_st114_clean <- snp_dist_clean_fn('snp_matrix_entcpx_st114.csv', 'CPO009')
snp_distances_entcpx_st190_clean <- snp_dist_clean_fn('snp_matrix_entcpx_st190.csv', 'CPO068')
snp_distances_entcpx_st93_clean <- snp_dist_clean_fn('snp_matrix_entcpx_st93.csv', 'CPO051')
snp_distances_koxy_st278_clean <- snp_dist_clean_fn('snp_matrix_koxy_st278.csv', 'CPO118')


snp_distances_all_clean <- bind_rows(snp_distances_serratia_lin_1_clean, snp_distances_serratia_lin_2_clean,
                                     snp_distances_pseudomonas_clean, snp_distances_entcpx_st114_clean,
                                     snp_distances_entcpx_st190_clean, snp_distances_entcpx_st93_clean,
                                     snp_distances_koxy_st278_clean)

snp_distances_all_clean_rev <- rename(snp_distances_all_clean, temp = isolate.x,
                                      isolate.x = isolate.y) %>%
  rename(isolate.y = temp)


#**All patient overlaps
pt_mvmt_isolate_metadata <- read_csv('pt_mvmt_isolate_metadata_20230613.csv')
cpe_pt_mvmt_3 <- read_csv('cpe_pt_mvmt_3_20230613.csv')

cpe_pt_mvmt_3_nodes <- distinct(cpe_pt_mvmt_3, ur_dummy)

cpe_pt_mvmt_isolates <- left_join(cpe_pt_mvmt_3, pt_mvmt_isolate_metadata, by = 'ur_dummy')
cpe_pt_mvmt_isolates_dummy <- cpe_pt_mvmt_isolates


#**Plasmid level
cpe_pt_mvmt_plasmid <- left_join(cpe_pt_mvmt_isolates, cpe_pt_mvmt_isolates_dummy, by = c("ward_code", 'rep_type_clean')) %>%
  filter(ur_dummy.x!=ur_dummy.y) %>%
  filter(!(rep_type_clean=='Chromosome'&species_st.x!=species_st.y)) %>%
  mutate(interval_1 = interval(transfer_in_ward_bed_date.x, transfer_out_ward_bed_date.x),
         interval_2 = interval(transfer_in_ward_bed_date.y, transfer_out_ward_bed_date.y),
         overlap = day(as.period(intersect(interval_1, interval_2), "days"))) %>%
  filter(!is.na(overlap)) %>%
  left_join(., snp_distances_all_clean, by = c("isolate.x", 'isolate.y')) %>%
  left_join(., snp_distances_all_clean_rev, by = c("isolate.x", 'isolate.y')) %>%
  mutate(transmission_type = ifelse(species_st.x == species_st.y, "Clonal", "Plasmid"),
         snp_dist = ifelse(is.na(snp_dist.x), snp_dist.y, snp_dist.x),
         snp_dist_10 = ifelse(snp_dist<11, "\u2264 10 SNVs", ">10 SNVs"),
         snp_dist_10 = ifelse(is.na(snp_dist_10), "\u2264 10 SNVs", snp_dist_10))


cpe_pt_mvmt_plasmid_2 <- select(cpe_pt_mvmt_plasmid, ur_dummy.x, ur_dummy.y, ward_code, rep_type_clean, species_st.x, species_st.y,
                                transfer_in_ward_bed_date.x, transmission_type, snp_dist_10, snp_dist)

cpe_pt_mvmt_plasmid_3 <- non_dup_fn_pt_mvmt(cpe_pt_mvmt_plasmid_2)

cpe_pt_mvmt_plasmid_3_edges <- select(cpe_pt_mvmt_plasmid_3, ur_dummy.x, ur_dummy.y, ward_code, transmission_type, snp_dist_10) %>%
  group_by(ward_code) %>%
  mutate(total_links = n(),
         ward_code_clean = ifelse(total_links>2, ward_code, 'Other'),
         ward_code_clean = ifelse(ward_code_clean=='4WA', 'Internal Medicine', ward_code_clean),
         ward_code_clean = ifelse(ward_code_clean=='5EA', 'Pulmonology', ward_code_clean),
         ward_code_clean = ifelse(ward_code_clean=='6EA', 'Burns', ward_code_clean),
         ward_code_clean = ifelse(ward_code_clean=='6WS', 'Surgery', ward_code_clean),
         ward_code_clean = ifelse(ward_code_clean=='7EA', 'Hematology', ward_code_clean),
         ward_code_clean = ifelse(ward_code_clean=='RC', 'Rehabilitation', ward_code_clean))

#Added this step to be able to facet into two groups
cpe_pt_mvmt_plasmid_4_edges <- filter(cpe_pt_mvmt_plasmid_3_edges, transmission_type == 'Clonal') %>%
  mutate(facet_type = "Clonal") %>%
  bind_rows(., cpe_pt_mvmt_plasmid_3_edges) %>%
  mutate(facet_type = ifelse(is.na(facet_type), "Clonal and plasmid", facet_type))

cpe_pt_mvmt_plasmid_5_edges <- filter(cpe_pt_mvmt_plasmid_3_edges, transmission_type == 'Plasmid') %>%
  mutate(facet_type = "Plasmid") %>%
  bind_rows(., cpe_pt_mvmt_plasmid_4_edges)

#Draw graph
cpe_pt_mvmt_network_plasmid <- graph_from_data_frame(cpe_pt_mvmt_plasmid_4_edges, directed=FALSE, vertices = cpe_pt_mvmt_3_nodes)
cpe_pt_mvmt_network_plasmid_2 <- ggraph(cpe_pt_mvmt_network_plasmid, layout = 'kk')+# layout = 'fr', weights = (edges_ben$Weight*1)) +
  geom_edge_link(aes(colour = factor(ward_code_clean),
                     width = cpe_pt_mvmt_plasmid_4_edges$transmission_type,
                     alpha = cpe_pt_mvmt_plasmid_4_edges$snp_dist_10)) +
  scale_edge_alpha_manual("SNV distance", values = c("\u2264 10 SNVs" = 1, ">10 SNVs" = 0.2)) +
  scale_edge_colour_brewer(name = 'Spatiotemporal overlap on ward', type = "seq", palette = "Set3", direction = 1) +
  scale_edge_width_manual("Transmission type", values = c("Clonal" = 1, "Plasmid" = 2)) +
  geom_node_point() +
  theme(strip.text = element_text(size = 14, face = 'bold'),
        legend.title = element_text(size = 14, face = 'bold'),
        legend.text = element_text(size= 12)) +
  facet_edges(~facet_type)

ggsave(plot = cpe_pt_mvmt_network_plasmid_2, 'cpe_pt_mvmt_network_plasmid_2_20220408.svg',
       width = 5500, height = 2500, unit = 'px')

#*Supp figure 1####
mobtyper_results_imp4_contigs_4_qtr <- metadata_final %>%
                                       arrange(acc_coll_date) %>% # arrange by date
                                       group_by(Year_quarter = quarter(acc_coll_date, with_year = FALSE)) %>% # Group by new quarter col by date
                                       ungroup() %>%
                                       mutate(YQ = as.yearqtr(paste(year,' Q',Year_quarter))) # Convert to graphable format


ridge_timeline_count_min <-   ggplot(data = mobtyper_results_imp4_contigs_4_qtr,
         aes(x = YQ,
             y = rep_type_clean,
             fill = species_st_final,
             colour = species_st_final,
             height = ..count..)) + # ..count.. works well if I want unscaled, generally useful for unfaceted plots. However for facets, use ..ndensity.. as this scales to 1.0
  geom_density_ridges(aes(color = species_st_final),
                      size = 0.5,
                      scale = 0.9,
                      panel_scaling = FALSE, # This means panel scaling consistent between facets..sooo useful
                      stat = "binline", binwidth = 0.25) + # can use binline or "identity" or "density"
  scale_fill_brewer(palette= "Paired") +
  scale_color_brewer(palette= "Paired") +
  scale_x_yearqtr(breaks = seq(2002.25, 2021, by=0.5), # quarters are literally considered decimals in R. Q1 = 0. Q2 = 0.25. Q4 = 0.75
                  limits = c(2002.25, 2021),
                  format = "%Y Q%q") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 14, face = 'bold'),
        axis.title.y = element_text(size = 14, face = 'bold'),
        strip.background = element_blank(),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12, face = 'bold'),
        legend.position = 'right',
        strip.text = element_blank()) +
  scale_y_discrete(position='right') +
  ylab("Genetic setting of blaIMP-4") +
  xlab("Year / Quarter") +
  facet_grid(rep_type_clean ~ .,
             scales = "free_y", space = "free", drop = TRUE) # This plots only lines with data without messing up the scale between facets. For this to work, panel_scaling = FALSE in geom_density_ridges

# Produce total histogram plot
g2 <- ggplot(mobtyper_results_imp4_contigs_4_qtr) +
  geom_histogram(aes(x = YQ), binwidth = 0.25, color = "black", fill = "light grey") +
  scale_x_yearqtr(breaks = seq(2002.25, 2021, by=0.5), # quarters are literally considered decimals in R. Q1 = 0. Q2 = 0.25. Q4 = 0.75
                  limits = c(2002.25, 2021),
                  format = "%Y Q%q") +
  ylab("Total no. of genomes") +
  scale_y_continuous(position='right') +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 14, face = 'bold'),
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank()
  )

ridge_timeline_count_min_histo <- ridge_timeline_count_min %>% insert_top(g2, height = 0.25)
ggsave(file="ridge_timeline_count_min_histo_20220510_revision.svg", plot=ridge_timeline_count_min_histo, width=4500, height=3750, units = "px")

#*Supp Figure 2A####
alfred_pseudomonas_tree_big <- alfred_tree_fn('pseudomonas_20221107.clean.core.aln.treefile', "", "CPO111")
ggsave(plot = alfred_pseudomonas_tree_big, file = 'alfred_pseudomonas_tree_big_20221214.svg', width = 1000, height = 2000, units = 'px')

#*Supp Figure 2A inset####
pseudomonas_tree <- chromo_tree_fn('pseudomonas_20221107.clean.core.aln.treefile', 1, expression(italic("")))
ggsave(plot = pseudomonas_tree, file = 'alfred_pseudomonas_tree_20221214.svg', height = 1000, width = 500, unit = 'px')

#*Supp Figure 2B####
alfred_st190_tree_big <- alfred_tree_fn('entcpx190.clean.core.aln.treefile', "CPO413-ENTCPX-02-21", "CPO068")
ggsave(plot = alfred_st190_tree_big, file = 'alfred_st190_tree_big_20220524.svg', width = 1000, height = 2000, units = 'px')

#*Supp Figure 2B inset####
alfred_st190_tree <- chromo_tree_fn('entcpx190.clean.core.aln.treefile', 1, expression(paste(italic(""), "")))
ggsave(plot = alfred_st190_tree, file = 'alfred_st190_tree_20220524.svg', height = 1000, width = 500, unit = 'px')

#*Supp Figure 2C####
alfred_st93_tree_big <- alfred_tree_fn('entcpx93.clean.core.aln.treefile', "CPO398-Entcpx-02-21", "CPO051")
ggsave(plot = alfred_st93_tree_big, file = 'alfred_st93_tree_big_20220524.svg', width = 1000, height = 2000, units = 'px')

#*Supp Figure 2C inset####
alfred_st93_tree <- chromo_tree_fn_2('entcpx93.clean.core.aln.treefile', 1, expression(paste(italic(""), ""))) #Toggled legend so that it would be present
ggsave(plot = alfred_st93_tree, file = 'alfred_st93_tree_20220524.svg', height = 1000, width = 700, unit = 'px')

#*Supp Figure 3####
flanker_figure_3_clusters <- read_csv('flanker_figure_3_clusters.csv') %>%
                              mutate(flanker_cluster = str_split_fixed(flanker_window_cluster,
                                                                       "_", n = 2)[, 2],
                                     flanker_cluster = as.double(flanker_cluster))

flanker_bubble_df <- select(metadata_final, isolate, rep_type_clean, flanker_cluster) %>%
                    filter(!is.na(flanker_cluster)) %>%
                    group_by(flanker_cluster) %>%
                    mutate(total_contigs = n()) %>%
                    filter(total_contigs>1) %>%
                    select(-total_contigs) %>%
                    group_by(rep_type_clean, flanker_cluster) %>%
                    summarise(total_flanker = n()) %>%
                    left_join(., flanker_figure_3_clusters, by = 'flanker_cluster') %>%
                    mutate(figure_3_cluster = ifelse(flanker_cluster==7, "H", figure_3_cluster),
                           figure_3_cluster = ifelse(flanker_cluster==8, "I", figure_3_cluster),
                           figure_3_cluster = factor(figure_3_cluster),
                           rep_type_clean = as.factor(rep_type_clean))

integron_snv_df <- select(metadata_final, isolate, rep_type_clean, nt_seq) %>%
                    filter(!is.na(nt_seq)) %>%
                    group_by(nt_seq) %>%
                    mutate(total_contigs = n()) %>%
                    filter(total_contigs>1) %>%
                    select(-total_contigs) %>%
                    group_by(rep_type_clean, nt_seq) %>%
                    summarise(total_nt_seq = n()) %>%
                    mutate(rep_type_clean = as.factor(rep_type_clean),
                           nt_seq = as.factor(nt_seq))

flanker_bubble_df_2 <- rename(flanker_bubble_df,
                              feature_var = figure_3_cluster, total_count = total_flanker) %>%
                        mutate(feature_var = paste('Cluster', feature_var, sep = " "),
                               facet_var = 'Flanking cluster')

integron_snv_df_2 <- rename(integron_snv_df,
                            feature_var = nt_seq, total_count = total_nt_seq) %>%
                      mutate(facet_var = 'Integron SNV group')

flanker_nt_seq_comb_df <- bind_rows(flanker_bubble_df_2,
                                    integron_snv_df_2)

flanker_nt_seq_bubble_plot <- ggplot(flanker_nt_seq_comb_df, aes(x=feature_var, y=rep_type_clean, size=total_count, colour = rep_type_clean)) + 
                              geom_point() +
                              geom_vline(xintercept = c('Cluster A', 'Cluster B', 'Cluster C', 'Cluster F',
                                                        'TGGTCGACGCCT' ,'GGGTCGACGCCT', 'GGGTCGACGTCT', 'GGGTCGACGTCT', 'GGGATGACGTCT'),
                                         linetype="dotted", size=0.5) +
                              scale_color_manual("Genetic setting", values = setNames(as.character(fill_source$fill_colour),
                                                                                      fill_source$`Genetic setting`), drop = TRUE) +                    
                              theme_classic() +
                              theme(axis.title.y = element_blank()) +
                              labs(size = "No. of genomes", color = 'Genetic setting') +
                              #facet_grid(cols = vars(facet_var), scales = 'free', space = 'free') #+
                              facet_grid(rows = vars(facet_var), scales = 'free', space = 'free') +
                              coord_flip() +
                              ylab("Genetic setting")

ggsave(plot = flanker_nt_seq_bubble_plot, file = 'flanker_nt_seq_bubble_plot_20221220.svg', width = 3200, height = 2600, unit = 'px')

