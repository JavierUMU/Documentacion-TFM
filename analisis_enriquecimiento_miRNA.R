#!/usr/bin/env Rscript

# Cargar librerías necesarias
library(dplyr)
library(readr)
library(GenomicRanges)
library(rtracklayer)
library(stringr)
library(gprofiler2)
library(ggplot2)
library(scales)

# Leer interacciones de miRTarBase y normalizar nombres
mirbase <- read_csv("miRTarBase_MTI.csv", col_types = cols()) %>%
  mutate(
    miRNA_low          = tolower(miRNA),
    species_miRNA_low  = tolower(`Species (miRNA)`),
    species_target_low = tolower(`Species (Target Gene)`)
  )

# Leer lista de variantes significativas
ids_significativos <- read_lines("variant_list.txt")

# Cargar datos de variantes y quedarnos con las que tienen coordenadas
df_variants <- readRDS("df_variants.rds")
df_sig <- df_variants %>%
  filter(
    id        %in% ids_significativos,
    !is.na(chr),
    !is.na(pos_start),
    !is.na(pos_end)
  )

# Crear objeto GRanges para las variantes
gr_variants <- GRanges(
  seqnames = df_sig$chr,
  ranges   = IRanges(start = df_sig$pos_start, end = df_sig$pos_end),
  id       = df_sig$id
)

# Importar regiones de miRNAs y buscar solapamientos
gr_mirnas <- import("mirnas.bed")
mirna_names <- mcols(gr_mirnas)$name
hits <- findOverlaps(gr_variants, gr_mirnas)
raw_hits <- unique(mirna_names[subjectHits(hits)])

# Extraer nombres de miRNAs legítimos
all_miRNAs <- str_split(raw_hits, "_") %>%
  unlist() %>%
  keep(~ str_detect(.x, regex("^hsa-miR-", ignore_case = TRUE))) %>%
  unique()
all_miRNAs_low <- tolower(all_miRNAs)

# Construir patrón para filtrar interacciones en miRTarBase
pattern_all <- paste0("^(", paste(all_miRNAs_low, collapse = "|"), ")(-3p|-5p)?$")

# Filtrar interacciones funcionales miRNA→gen
genes_diana <- mirbase %>%
  filter(
    str_detect(miRNA_low, pattern_all),
    species_miRNA_low  == "hsa",
    species_target_low == "hsa",
    `Support Type`     == "Functional MTI"
  ) %>%
  select(miRNA, `Target Gene`, `Target Gene (Entrez ID)`, Experiments) %>%
  distinct()

# Resumir número y lista de genes por miRNA base
resumen_por_miRNA <- genes_diana %>%
  mutate(miRNA_base = str_remove(tolower(miRNA), "-(3p|5p)$")) %>%
  group_by(miRNA_base) %>%
  summarise(
    n_genes = n_distinct(`Target Gene`),
    genes   = paste(unique(`Target Gene`), collapse = "; "),
    .groups = "drop"
  )

# Preparar vector de IDs de genes para g:Profiler
genes_unicos_entrez <- genes_diana %>%
  pull(`Target Gene (Entrez ID)`) %>%
  na.omit() %>%
  unique() %>%
  as.character()

# Ejecutar enriquecimiento con g:Profiler
gostres <- gost(
  query             = genes_unicos_entrez,
  organism          = "hsapiens",
  correction_method = "fdr",
  sources           = c("GO:BP", "GO:MF", "KEGG", "REAC", "HP"),
  user_threshold    = 0.05,
  evcodes           = TRUE
)

# Extraer resultados y calcular -log10(p-value)
gost_df <- gostres$result %>%
  mutate(logp = -log10(p_value))

# Filtrar términos HPO relacionados con corazón
heart_hp <- gost_df %>%
  filter(source == "HP", str_detect(term_name, regex("card|heart", ignore_case = TRUE))) %>%
  arrange(intersection_size) %>%
  mutate(term_name = factor(term_name, levels = term_name))

# Contar términos totales y GO:BP de corazón
total_terms  <- nrow(gost_df)
heart_bp     <- gost_df %>%
  filter(source == "GO:BP", str_detect(term_name, regex("card|heart", ignore_case = TRUE))) %>%
  nrow()

cat("Número total de términos:", total_terms, "\n")
cat("Términos GO:BP relacionados con corazón:", heart_bp, "\n")

# Filtrar top 30 términos KEGG
top30_kegg <- gost_df %>%
  filter(source == "KEGG") %>%
  slice_min(order_by = p_value, n = 30) %>%
  arrange(p_value) %>%
  mutate(term_name = factor(term_name, levels = term_name))

# Gráfico de burbujas para KEGG
ggplot(top30_kegg, aes(x = logp, y = term_name, size = intersection_size)) +
  geom_point(alpha = 0.7) +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed") +
  scale_size_area(max_size = 8) +
  labs(
    x     = expression(-log[10](p-value)),
    y     = NULL,
    size  = "Genes en ruta",
    title = "Top 30 términos KEGG"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.y  = element_text(size = 6),
    plot.title   = element_text(size = 12, hjust = 0.5)
  )

# Filtrar top 30 términos HPO
top30_hp <- gost_df %>%
  filter(source == "HP") %>%
  slice_min(order_by = p_value, n = 30) %>%
  arrange(p_value) %>%
  mutate(term_name = factor(term_name, levels = term_name))

# Gráfico de burbujas para HPO
ggplot(top30_hp, aes(x = logp, y = term_name, size = intersection_size)) +
  geom_point(alpha = 0.7, color = "#e31a1c") +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed") +
  scale_size_area(max_size = 8) +
  labs(
    x     = expression(-log[10](p-value)),
    y     = NULL,
    size  = "Genes asociados",
    title = "Top 30 términos HPO"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.y  = element_text(size = 6),
    plot.title   = element_text(size = 12, hjust = 0.5)
  )

# Gráfico de barras para términos HPO de corazón
ggplot(heart_hp, aes(x = intersection_size, y = term_name, fill = p_value)) +
  geom_col(color = "black", width = 0.7) +
  scale_fill_gradient(low = "firebrick", high = "lightpink", name = "p-value") +
  labs(
    x     = "Número de genes",
    y     = NULL,
    title = "Términos HPO relacionados con cardiopatías"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.y = element_text(size = 8, face = "bold"),
    plot.title  = element_text(size = 14, hjust = 0.5)
  ) +
  coord_cartesian(expand = FALSE)
