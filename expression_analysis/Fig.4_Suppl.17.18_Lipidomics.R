library(lipidr)
library(dplyr)
library(limma)
library(rgoslin)
library(ggplot2)
`%not in%` <- Negate(`%in%`)

# Data ----
dir <- "Lipidomics"

meta <- read.csv(file.path(dir, "Lipidomics_Meta.tsv"),
                 sep = "\t")
df <- read.csv(file.path(dir, "Lipidomics_Data.tsv"),
               sep = "\t")
df_params <- read.csv(file.path(dir, "Lipidomics_LipidParameters.tsv"),
                      sep = "\t")
df_params$Order <- seq(9, nrow(df_params) + 8, 1)
df$Order <- seq(9, nrow(df) + 8, 1)
df_merge <- merge(df, df_params %>% select(Order, Lipid), by = "Order", all.y = T)

# add lipid names
df_merge$Lipid_both <- df_merge$Lipid
df_merge$Lipid <- gsub("\\|.*", "", df_merge$Lipid)
df_merge$Lipid <- gsub("_", "/", df_merge$Lipid)
df_merge$Lipid <- sub(" O-", "O ", df_merge$Lipid)
df_merge$Lipid <- sub(" P-", "P ", df_merge$Lipid)
df_merge$Lipid <- gsub(";(O\\d*)", "(\\1)", df_merge$Lipid)
df_merge$Lipid <- gsub(";(\\d*O)", "(\\1)", df_merge$Lipid)
df_merge$Lipid <- gsub(";(S)", "(\\1)", df_merge$Lipid)
df_merge <- df_merge %>%
  select(Lipid, lipids, colnames(.)[3:26])
not_uniq <- (df_merge %>% group_by(Lipid) %>%
               mutate(lip_n = rep_len(1:10, length.out = n())) %>%
               summarise(lip_n_sum = sum(lip_n)) %>%
               filter(lip_n_sum > 1))$Lipid
df_merge <- df_merge %>%
  group_by(Lipid) %>%
  mutate(lip_n = rep_len(1:10, length.out = n()))
df_merge$Lipid <- ifelse(df_merge$Lipid %in% not_uniq,
                         sprintf("%s (%s)", df_merge$Lipid, df_merge$lip_n),
                         df_merge$Lipid)

# Analysis (by condition) ----
## separately for 'all'/'FrL'/'Hyp'

reg <- "all"

d <- as_lipidomics_experiment(df_merge %>% select(-c("lipids", "lip_n")))
d <- add_sample_annotation(d, annot_file = meta)
coq <- c('CoQ9', 'CoQ10')
rowData(d[coq, ])$Class <- 'CoQ'

if (reg != "all") {d = d[, d$Region == reg]}

## PCA ----
# Supplementary Figure 17A
mvaresults = mva(d, measure="Area", method="PCA")
df <- data.frame(PC1 = mvaresults$scores$p1,
                 PC2 = mvaresults$scores$p2)
rownames(df) <- mvaresults[["col_data"]]@rownames
df <- merge(df %>% tibble::rownames_to_column("Sample"),
            meta, by = "Sample")
ggplot(df, aes(x = PC1, y = PC2, fill = Group)) +
  geom_point(aes(shape = Sex), size = 4, stroke = 0.5) +  
  scale_shape_manual(values = c(21, 24)) +  
  scale_fill_manual(
    values = c("#56B1F7", "#fb7669"),
    labels = c("Ctrl", "5HTP")
  ) +
  theme_minimal() +
  theme(legend.position = "right", 
        axis.line = element_line(size = 0.15),
        axis.ticks = element_line(size = 0.1),
        legend.title = element_text(size = 18, color = "black"),
        legend.text = element_text(size = 15, color = "black"),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 15, color = "black")) +
  guides(
    shape = guide_legend(
      override.aes = list(fill = "white", color = "black") 
    ),
    fill = guide_legend(
      override.aes = list(shape = 21, color = "black")
    )
  )


## DE ----
de_results <- de_analysis(
  data=d, group_col = "Group",
  HT - Ctrl,
  measure="Area"
)

## Enrichment ----
de_results$Class <- gsub("MGDGO", "MGDG", de_results$Class)
enrich_results = lsea(de_results, rank.by = "logFC", min_size = 5)

enrich_sign <- enrich_results %>% filter(padj < 0.05) %>% 
  mutate(Ctrl = ifelse(NES < 0, T, F),
         HT = ifelse(NES > 0, T, F)) %>% 
  select(set, Ctrl, HT) %>% 
  mutate(class = gsub("Class_", "", set))
df <- de_results %>% 
  select(Class, logFC) %>% 
  mutate(Significant = case_when(
    Class %in% enrich_sign$class[enrich_sign$Ctrl] ~ "Ctrl",
    Class %in% enrich_sign$class[enrich_sign$HT] ~ "5HTP",
    TRUE ~ "not_sign"
  ))

# Figure 4J / Supplementary Figure 18A
df$Class <- factor(df$Class,
                   levels = rev(c('PC', 'PE', 'LPC', 'LPE', 'PI', 'PS', 'PG', 
                                  'BMP', "HBMP", 'CL', 'TG', 'DG', 'MGDG', 'HexCer', 
                                  'Cer', 'SM', 'FA', 'CAR', 'ST', 'CoQ')))
ggplot(df, aes(y = Class, x = logFC, fill = Significant)) +
  geom_boxplot(color = "grey40", lwd=0.3, outlier.size = 1) + 
  geom_vline(xintercept = 0, lty = 2) +
  scale_fill_manual(
    values = c(`not_sign` = "white", `Ctrl`="#56B1F7", `5HTP` = "#fb7669"), 
    drop=FALSE) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.line = element_line(size = 0.2),
        axis.text.x = element_text(size = 12, color = "black", angle = 90),
        axis.text.y = element_text(size = 14, color = "black"),
        axis.title.y = element_blank(),
        axis.ticks.x = element_line(),
        axis.title.x = element_text(size = 15)) 

# Figure 4L / Supplementary Figure 18D
plot_list <- list()
for (pthws in c("MGDG", "PC", "PE")) {
  lipids <- (de_results %>% 
               filter(Molecule %in% grep(paste0("^", pthws), de_results$Molecule, value = T))
  )$Molecule
  
  p <- fgsea::plotEnrichment(lipids,
                        setNames(de_results$logFC, de_results$Molecule), 
                        ticksSize = 0.7) +
    labs(title=pthws) +
    xlab("rank of lipids") +
    theme(plot.margin = margin(30,5,5,5))
  plot_list[[pthws]] <- p

}
cowplot::plot_grid(plotlist = plot_list, ncol = 1)


## Intensities ----
### Sphingomyelins ----
# Figure 4M / Supplementary Figure 18B
sm <- grep("SM ", de_results$Molecule, value = T); sm
sm_sign <- grep("SM ", (de_results %>% filter(P.Value < 0.05))$Molecule, value = T); sm_sign
df <- d@assays@data@listData[["Area"]] %>% as.data.frame() %>% 
  tibble::rownames_to_column("lipid") %>% 
  filter(lipid %in% sm)
df_plot <- df %>% tidyr::gather(key = "sample", value = "Intensity", -lipid) %>% 
  mutate(group = gsub(".*(Ctrl|HT).*", "\\1", sample))
df_plot$group <- factor(gsub("HT", "5HTP", df_plot$group),
                        levels = c("Ctrl", "5HTP"))
df_plot$lipid <- ifelse(df_plot$lipid %in% sm_sign, 
                        paste0(df_plot$lipid, " *"),
                        df_plot$lipid)

ggplot(df_plot, aes(x = group, y = Intensity, fill = group)) +
  geom_boxplot() +
  geom_point() +
  facet_wrap(lipid ~ ., scales = "free", nrow = 2) +
  scale_fill_manual(values = c("#56B1F7", "#fb7669")) +
  theme(panel.background = element_blank(),
        panel.border = element_rect(linewidth = 0.5, color = "grey80", fill = NA),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 9, face = "bold"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 18, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"))

### Cholesterol ----
# Figure 4I / Supplementary Figure 18C
df <- d@assays@data@listData[["Area"]] %>% as.data.frame() %>% 
  tibble::rownames_to_column("lipid") %>% 
  filter(lipid %in% c("ST 27:1(O) (Cholesterol)"))
df_plot <- df %>% tidyr::gather(key = "sample", value = "Intensity", -lipid) %>% 
  mutate(group = gsub(".*(Ctrl|HT).*", "\\1", sample))
df_plot$lipid <- gsub(" (Cholesterol)", "\n(Cholesterol)", df_plot$lipid)
df_plot$group <- factor(gsub("HT", "5HTP", df_plot$group),
                        levels = c("Ctrl", "5HTP"))
ggplot(df_plot %>% 
         filter(grepl("Cholesterol", lipid)), 
       aes(x = group, y = Intensity, fill = group)) +
  geom_boxplot() +
  geom_point() +
  facet_wrap(lipid ~ ., scales = "free") +
  scale_fill_manual(values = c("#56B1F7", "#fb7669")) +
  theme(panel.background = element_blank(),
        panel.border = element_rect(linewidth = 0.5, color = "grey80", fill = NA),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 12, face = "bold"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 15, color = "black"),
        axis.title.y = element_text(size = 18, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"))




# Analysis (by region) ----
## separately for 'Ctrl' / 'HT' / 'all'

cond <- "HT"

d <- as_lipidomics_experiment(df_merge %>% select(-c("lipids", "lip_n")))
d <- add_sample_annotation(d, annot_file = meta)
coq <- c('CoQ9', 'CoQ10')
rowData(d[coq, ])$Class <- 'CoQ'

if (cond != "all") {d = d[, d$Group == cond]}


## PCA ----
# Supplementary Figure 17C/E

if (cond == "all") {
  print("see PCA from the previous section")
} else {
  mvaresults = mva(d, measure="Area", method="PCA")
  df <- data.frame(PC1 = mvaresults$scores$p1,
                   PC2 = mvaresults$scores$p2)
  rownames(df) <- mvaresults[["col_data"]]@rownames
  df <- merge(df %>% tibble::rownames_to_column("Sample"),
              meta, by = "Sample")
  ggplot(df, aes(x = PC1, y = PC2, fill = Region)) +
    geom_point(aes(shape = Sex), size = 4, stroke = 0.5) +  
    scale_shape_manual(values = c(21, 24)) +  
    scale_fill_manual(
      values = c("goldenrod", "tomato3"),
      labels = c("Cortex", "Hypothalamus")
    ) +
    theme_minimal() +
    theme(legend.position = "right", 
          axis.line = element_line(size = 0.15),
          axis.ticks = element_line(size = 0.1),
          legend.title = element_text(size = 18, color = "black"),
          legend.text = element_text(size = 15, color = "black"),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 15, color = "black")) +
    guides(
      shape = guide_legend(
        override.aes = list(fill = "white", color = "black")
      ),
      fill = guide_legend(
        override.aes = list(shape = 21, color = "black")  
      )
    )
}


## DE ----
de_results <- de_analysis(
  data=d, group_col = "Region",
  Hyp - FrL,
  measure="Area"
)
# Supplementary Figure 18B/D/F
EnhancedVolcano::EnhancedVolcano(de_results, 
                                 pCutoff = 0.05, FCcutoff = 0.1,
                                 lab = de_results$Molecule,
                                 x = 'logFC',
                                 y = 'adj.P.Val',
                                 title = reg,
                                 subtitle = "Cortex - Hypothalamus")

## Enrichment ----
de_results$Class <- gsub("MGDGO", "MGDG", de_results$Class)
enrich_results = lsea(de_results, rank.by = "logFC", min_size = 5)

enrich_sign <- enrich_results %>% filter(padj < 0.05) %>% 
  mutate(Cortex = ifelse(NES < 0, T, F),
         Hyp = ifelse(NES > 0, T, F)) %>% 
  select(set, Cortex, Hyp) %>% 
  mutate(class = gsub("Class_", "", set))
df <- de_results %>% 
  select(Class, logFC) %>% 
  mutate(Significant = case_when(
    Class %in% enrich_sign$class[enrich_sign$Cortex] ~ "Cortex",
    Class %in% enrich_sign$class[enrich_sign$Hyp] ~ "Hyp",
    TRUE ~ "not_sign"
  ))

# Supplementary Figure 18I
df$Class <- factor(df$Class,
                   levels = rev(c('PC', 'PE', 'LPC', 'LPE', 'PI', 'PS', 'PG', 
                                  'BMP', "HBMP", 'CL', 'TG', 'DG', 'MGDG', 'HexCer', 
                                  'Cer', 'SM', 'FA', 'CAR', 'ST', 'CoQ')))

ggplot(df, aes(y = Class, x = logFC, fill = Significant)) +
  geom_boxplot(color = "grey40", lwd=0.3, outlier.size = 1) + 
  geom_vline(xintercept = 0, lty = 2) +
  scale_fill_manual(
    values = c(`not_sign` = "white", `Cortex`="goldenrod", `Hyp` = "tomato3"), 
    drop=FALSE) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.line = element_line(size = 0.2),
        axis.text.x = element_text(size = 12, color = "black", angle = 90),
        axis.text.y = element_text(size = 14, color = "black"),
        axis.title.y = element_blank(),
        axis.ticks.x = element_line(),
        axis.title.x = element_text(size = 15)) 

# Supplementary Figure 18J
plot_list <- list()
for (pthws in intersect(enrich_sign$class, df_params$Class)) {
  lipids <- (de_results %>% 
               filter(Molecule %in% grep(paste0("^", pthws), de_results$Molecule, value = T))
  )$Molecule
  
  p <- fgsea::plotEnrichment(lipids,
                        setNames(de_results$logFC, de_results$Molecule), 
                        ticksSize = 0.7) +
    labs(title=pthws) +
    xlab("rank of lipids") +
    theme(plot.margin = margin(20,5,5,5))
  plot_list[[pthws]] <- p
}
cowplot::plot_grid(plotlist = plot_list, nrow = 1)


## Compare Ctrl - HT ----

# calculted in previous sections
de_ctrl <- read.csv(file = "ctrl_de_lipidr.csv") %>% 
  filter(adj.P.Val < 0.05)
de_ht <- read.csv(file = "ht_de_lipidr.csv") %>% 
  filter(adj.P.Val < 0.05)

### Cortex (down)
# Supplementary Figure 18G (left)
lipid_list <- list(Ctrl = (de_ctrl %>% filter(logFC < 0))$Molecule,
                   HTP = (de_ht %>% filter(logFC < 0))$Molecule)
names(lipid_list) <- c("Ctrl", "5HTP")
ggVennDiagram::ggVennDiagram(lipid_list, label = "count", 
              set_size = 3, label_size = 8) + 
  scale_fill_distiller(palette = "Reds", direction = 1) +
  scale_x_continuous(expand = expansion(mult = 0.2)) +
  theme(legend.key.size = unit(0.8, "cm"),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 15)) +
  ggtitle("Cortex") +
  theme(plot.title = element_text(size = 22, hjust = 0.5))


### Hyp (up)
# Supplementary Figure 18G (right)
lipid_list <- list(Ctrl = (de_ctrl %>% filter(logFC > 0))$Molecule,
                   HTP = (de_ht %>% filter(logFC > 0))$Molecule)
names(lipid_list) <- c("Ctrl", "5HTP")
ggVennDiagram::ggVennDiagram(lipid_list, label = "count", 
              set_size = 3, label_size = 8) + 
  scale_fill_distiller(palette = "Reds", direction = 1) +
  scale_x_continuous(expand = expansion(mult = 0.2)) +
  theme(legend.key.size = unit(0.8, "cm"),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 15)) +
  ggtitle("Hypothalamus") +
  theme(plot.title = element_text(size = 22, hjust = 0.5))



## Comparison of differentially abundant lipids in brain regions in different conditions
# Supplementary Figure 18H

### Cortex (down)
lipid_list <- list(Ctrl = (de_ctrl %>% filter(logFC < 0))$Class,
                   HTP = (de_ht %>% filter(logFC < 0))$Class)
names(lipid_list) <- c("Ctrl", "5HTP")

df1 <- data.frame(Region = "Cortex",
                  Condition = c(rep("Ctrl", length(lipid_list$Ctrl)),
                                rep("5HTP", length(lipid_list$`5HTP`))),
                  Lipids = c(lipid_list$Ctrl, lipid_list$`5HTP`)
)

### Hyp (up)
lipid_list <- list(Ctrl = (de_ctrl %>% filter(logFC > 0))$Class,
                   HTP = (de_ht %>% filter(logFC > 0))$Class)
names(lipid_list) <- c("Ctrl", "5HTP")

df2 <- data.frame(Region = "Hypothalamus",
                  Condition = c(rep("Ctrl", length(lipid_list$Ctrl)),
                                rep("5HTP", length(lipid_list$`5HTP`))),
                  Lipids = c(lipid_list$Ctrl, lipid_list$`5HTP`)
)

### all
df <- rbind(df1, df2)

df$Lipids <- factor(df$Lipids,
                    levels = c('PC', 'PE', 'LPC', 'LPE', 'PI', 'PS', 'PG', 
                               'BMP', "HBMP", 'CL', 'TG', 'DG', 'MGDG', 'HexCer', 
                               'Cer', 'SM', 'FA', 'CAR', 'ST', 'CoQ'))

ggplot(df, aes(x = Lipids, group = Condition, fill = Condition)) +
  geom_bar(position = position_fill()) +
  facet_grid(Region ~ .) +
  theme_minimal() +
  theme(axis.title = element_blank(),
        axis.text.x = element_text(size = 18, angle = 45, hjust = 1, color = "black"),
        axis.text.y = element_text(size = 8), 
        strip.text = element_text(size = 10, face = "bold"),
        panel.grid = element_blank()) +
  scale_fill_manual(values = rev(c("#56B1F7", "#fb7669")))