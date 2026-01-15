#LIBRARIES
library(ggplot2)
library(readr)
library(dplyr)
library(tidyr)
library(factoextra)
library(FactoMineR)
library(ggrepel)
library(pheatmap)
library(RColorBrewer)
library(patchwork)
library(ggbreak)
library(tibble)



# COLOR PALETTE
custom_colors <- c(
  "California"   = "#1f77b4",
  "Florida"      = "#ff7f0e",
  "Hawaii"       = "#9467bd",
  "New York"     = "#d62728",
  "Puerto Rico"  = "#2ca02c"
)
custom_colors_all <- c(custom_colors, "Total" = "goldenrod2")

# DIRECTORIES
base_dir  <- "/Users/path"
plots_dir <- file.path(base_dir, "plots")
if (!dir.exists(plots_dir)) dir.create(plots_dir)

# FIGURE 1 — BAYESIAN SKYLINE PLOT (EPS, 400 dpi)
skyline_paths <- list(
  California  = "bayesian_skyline_ca.txt",
  Florida     = "bayesian_skyline_fl.txt",
  Hawaii      = "bayesian_skyline_hi.txt",
  NewYork     = "bayesian_skyline_ny.txt",
  PuertoRico  = "bayesian_skyline_pr.txt"
)
sky_dir <- file.path(base_dir, "Tracer_log_file")

skyline_dfs <- lapply(names(skyline_paths), function(loc) {
  df <- read.table(file.path(sky_dir, skyline_paths[[loc]]), header = TRUE, sep = "", comment.char = "")
  if (ncol(df) >= 8) colnames(df)[1:8] <- c("time","date","datetime","ms","mean","median","upper","lower")
  df %>% select(time,mean,median,lower,upper) %>% mutate(location = loc)
})
combined_skyline <- bind_rows(skyline_dfs) %>%
  mutate(location = recode(location,
                           California="California", Florida="Florida",
                           Hawaii="Hawaii", NewYork="New York",
                           PuertoRico="Puerto Rico"))

p_skyline <- ggplot(combined_skyline,aes(x=time,y=median,color=location))+
  geom_line(linewidth=1.2)+
  geom_ribbon(aes(ymin=lower,ymax=upper,fill=location),alpha=0.25,color=NA)+
  scale_color_manual(values=custom_colors)+
  scale_fill_manual(values=custom_colors)+
  theme_classic(base_size=14)+
  labs(y=expression("Viral Effective Population Size (Ne)"),
       x="",color="Location",fill="Location")+
  geom_vline(xintercept=2020.25,color="gray40",linewidth=0.7)+
  geom_vline(xintercept=2020.45,color="gray40",linewidth=0.7)+
  geom_vline(xintercept=2021.0 ,color="gray40",linewidth=0.7)+
  annotate("text",x=2020.25,y=220,label="Stay-at-home declared",angle=90,vjust=-0.4,size=3.2)+
  annotate("text",x=2020.45,y=220,label="Stay-at-home order lifted",angle=90,vjust=-0.4,size=3.2)+
  annotate("text",x=2021.0 ,y=220,label="Vaccination starts",angle=90,vjust=-0.4,size=3.2)

ggsave(file.path(plots_dir,"Figure1_BayesianSkyline_SOLIDLINES.eps"),
       plot=p_skyline,width=10,height=7,dpi=400,device=cairo_ps)

#FIGURE 2 — Positive Selection (A) + Epistasis (B) (EPS, 400 dpi)
## A. Positive Selection 
sig_data <- read_csv(file.path(base_dir,"POSITIVE_selection_vision_HyPhy/significant_sites_table.csv"))
protein_map <- data.frame(
  Protein=c("nsp1","nsp2","nsp3","nsp4","nsp5","nsp6","nsp7","nsp8","nsp9","nsp10",
            "nsp11","nsp12","nsp13","nsp14","nsp15","nsp16","S","ORF3a","E","M",
            "ORF6","ORF7a","ORF7b","ORF8","N","ORF10"),
  Start=c(1,181,819,2764,3264,3569,3854,3943,4141,4254,
          4393,4394,5325,5926,6453,6797,7097,8555,8803,8943,
          10020,10152,10323,10438,10734,10901),
  End=c(180,818,2763,3263,3568,3853,3942,4140,4253,4392,
        4394,5324,5925,6452,6796,7096,8554,8802,8942,10019,
        10151,10322,10437,10733,10900,11055)
)
positive_sites <- sig_data %>%
  filter(`p-value`<0.05)%>%
  rowwise()%>%
  mutate(Protein=protein_map$Protein[which(Site>=protein_map$Start & Site<=protein_map$End)])%>%
  ungroup()
protein_location_counts <- positive_sites %>%
  group_by(Protein,Location)%>%
  summarise(n_sites=n(),.groups="drop")
totals_pos <- protein_location_counts %>%
  group_by(Protein)%>%
  summarise(value=sum(n_sites))%>%
  mutate(type="Total")
regions_pos <- protein_location_counts %>%
  mutate(value=-n_sites)%>%
  rename(type=Location)%>%
  select(Protein,value,type)
ordered_proteins_pos <- totals_pos %>% arrange(desc(value)) %>% pull(Protein)
totals_pos$Protein  <- factor(totals_pos$Protein ,levels=ordered_proteins_pos)
regions_pos$Protein <- factor(regions_pos$Protein,levels=ordered_proteins_pos)
legend_order <- c("California","Florida","New York","Hawaii","Puerto Rico","Total")
regions_pos$type <- factor(regions_pos$type,levels=legend_order)
totals_pos$type  <- factor(totals_pos$type ,levels=legend_order)

combi_pos <- ggplot()+
  geom_bar(data=totals_pos ,aes(x=Protein,y=value,fill=type),stat="identity",width=0.8)+
  geom_bar(data=regions_pos,aes(x=Protein,y=value,fill=type),
           stat="identity",position_dodge(width=0.9),width=0.7)+
  scale_fill_manual(values=custom_colors_all)+
  scale_y_continuous(breaks=pretty(c(-40,90)),labels=function(x)abs(x))+
  theme_classic(base_size=12)+
  labs(x="",y="Count of Positively Selected Sites",fill="Location")+
  geom_hline(yintercept=0,color="black")+
  theme(axis.text.x=element_text(angle=45,hjust=1),
        legend.position="top",
        legend.title=element_text(face="bold"),
        legend.text=element_text(size=10))+
  guides(fill=guide_legend(nrow=1))

## B. Epistasis 
epi_data <- read.csv(file.path(base_dir,"EIs_datasets/epistatic_sites_per_protein_per_region.csv"))
epi_data <- epi_data[!epi_data$protein %in% c("3prime","NegStrand"),]
totals_epi <- epi_data %>% group_by(protein)%>% summarise(value=sum(n))%>% mutate(type="Total")
regions_epi<- epi_data %>% mutate(value=-n)%>% rename(type=region)%>% select(protein,value,type)
ordered_proteins_epi <- totals_epi %>% arrange(desc(value))%>% pull(protein)
totals_epi$protein  <- factor(totals_epi$protein ,levels=ordered_proteins_epi)
regions_epi$protein <- factor(regions_epi$protein,levels=ordered_proteins_epi)
combi_p <- ggplot()+
  geom_bar(data=totals_epi ,aes(x=protein,y=value,fill=type),stat="identity",width=0.8)+
  geom_bar(data=regions_epi,aes(x=protein,y=value,fill=type),
           stat="identity",position_dodge(width=0.9),width=0.7)+
  scale_fill_manual(values=custom_colors_all)+
  scale_y_continuous(breaks=pretty(c(-700,2300)),labels=function(x)abs(x))+
  scale_y_break(c(500,2000),scales=0.5)+
  theme_classic(base_size=12)+
  labs(x="Protein",y="Number of Epistatic Sites")+
  geom_hline(yintercept=0,color="black")+
  theme(axis.text.x=element_text(angle=45,hjust=1),
        legend.position="none")

## C. Figura 2 combinada A + B 
final_plot <- (combi_pos / combi_p) +
  plot_layout(guides="collect") +
  plot_annotation(tag_levels="A") &
  theme(legend.position="bottom",
        plot.tag=element_text(face="bold",size=14))

ggsave(file.path(plots_dir,"Figure2_PositiveSelection_Epistasis_COMBINED.eps"),
       plot=final_plot,width=14,height=10,dpi=400,device=cairo_ps)

#FIGURE 3 — PCA (EPS, 400 dpi)


# Leer y limpiar datos
epistasis_df <- read_csv(file.path(base_dir, "EIs_datasets/epistasis_summary_for_pca.csv")) %>%
  filter(!is.na(region), region != "") %>%
  select(region, num_nodes, num_edges, avg_degree, density, num_components, high_degree_nodes)

# Preparar matriz para PCA
epistasis_data <- epistasis_df %>%
  column_to_rownames("region") %>%
  select(where(is.numeric))

res.pca <- PCA(epistasis_data, graph = FALSE)

# Extraer coordenadas
ind_df <- as.data.frame(res.pca$ind$coord)
ind_df$region <- rownames(ind_df)
region_labels <- c(
  ca = "California",
  fl = "Florida",
  hi = "Hawaii",
  ny = "New York",
  pr = "Puerto Rico"
)
ind_df$label <- region_labels[ind_df$region]
ind_df$label[is.na(ind_df$label)] <- ind_df$region  # por si hay algo extra
ind_df$color <- custom_colors[ind_df$label]

# PCA plot con colores y círculos sólidos
gg_pca <- ggplot(ind_df, aes(x = Dim.1, y = Dim.2)) +
  geom_point(aes(color = label), size = 5, shape = 21, stroke = 1, fill = NA) +
  geom_point(aes(fill = label), size = 5, shape = 21, color = "black") +
  scale_fill_manual(values = custom_colors) +
  scale_color_manual(values = custom_colors) +
  geom_text_repel(aes(label = label), size = 4.5, show.legend = FALSE) +
  theme_classic(base_size = 14) +
  labs(
    x = paste0("PC1 (", round(res.pca$eig[1, 2], 1), "%)"),
    y = paste0("PC2 (", round(res.pca$eig[2, 2], 1), "%)")
  ) +
  theme(legend.position = "none")

# Guardar figura en EPS (400 ppi)
ggsave(
  file.path(plots_dir, "Figure3_PCA_epistasis_summary.eps"),
  plot = gg_pca,
  width = 8,
  height = 6,
  dpi = 400,
  device = cairo_ps
)


# FIGURE 4 — Heatmap (TIFF, 500 dpi)


# Leer datos
heat_data <- read.csv(file.path(base_dir, "EIs_datasets/epistatic_sites_per_protein_per_region.csv"))

# Excluir elementos no deseados
heat_data <- heat_data %>%
  filter(!protein %in% c("3prime", "NegStrand"))

# Orden genómico (basado en SARS-CoV-2)
genomic_order <- c(
  "nsp1", "nsp2", "nsp3", "nsp4", "nsp5", "nsp6",
  "nsp7", "nsp8", "nsp9", "nsp10", "nsp11", "nsp12",
  "nsp13", "nsp14", "nsp15", "nsp16",
  "S", "ORF3a", "E", "M", "ORF6", "ORF7a",
  "ORF7b", "ORF8", "N", "ORF10"
)

# Reordenar y pivotear datos (regiones en columnas)
heat_matrix <- heat_data %>%
  mutate(protein = factor(protein, levels = genomic_order)) %>%
  arrange(protein) %>%
  pivot_wider(names_from = region, values_from = n, values_fill = 0) %>%
  column_to_rownames("protein")

# Calcular Z-score por fila (proteína)
z_mat <- t(scale(t(as.matrix(heat_matrix))))

# Asegurar orden deseado en columnas (regiones)
desired_order <- c("California", "Florida", "New York", "Puerto Rico", "Hawaii")
z_mat <- z_mat[, desired_order[desired_order %in% colnames(z_mat)]]

# Paleta personalizada: violeta ↔ blanco ↔ amarillo
heat_colors <- colorRampPalette(c("goldenrod2", "white", "#1f77b4"))(100)

# Exportar heatmap
tiff(file.path(plots_dir, "Figure4_Heatmap_epistasis_zscore_horizontal_labels.tiff"),
     width = 8, height = 7, units = "in", res = 500)

pheatmap(
  z_mat,
  color = heat_colors,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  border_color = NA,
  fontsize = 10,
  angle_col = 0,        
  main = ""
)

dev.off()
