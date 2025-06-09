library(allelic)
library(mongolite)
library(dplyr)
library(tidyr)
library(tibble)

uri <- "mongodb://localhost:27017/?"

# Leer colecciones desde MongoDB
df_variants <- mongo("variants",          "MIRNAS", url = uri)$find()
df_sv       <- mongo("samples_variants", "MIRNAS", url = uri)$find()
df_sample   <- mongo("samples",           "MIRNAS", url = uri)$find()
df_effects  <- mongo("effects",           "MIRNAS", url = uri)$find()

# Filtrar muestras y definir grupos de interés
canalo_list <- c(
  "Alt. Conducción", "Brugada", "QT largo",
  "SADS posible (SUDS)", "SADS-corazon normal", "TV catecol"
)

df_s <- df_sample %>%
  filter(patol_prin %in% c("Hipertrofica", canalo_list)) %>%
  transmute(
    sample    = id,
    status    = if_else(patol_prin == "Hipertrofica", "case", "control"),
    patologia = if_else(status == "case", "Hipertrofia", "Canalopatias")
  )

# Unir genotipos con estado de la muestra
df_join <- df_sv %>%
  inner_join(df_s, by = "sample")

# Calcular conteos homocigoto y heterocigoto por variante y grupo
counts2 <- df_join %>%
  group_by(variant, status) %>%
  summarise(
    hom_alt = sum(homo, na.rm = TRUE),
    het     = sum(!homo, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(
    df_s %>% count(status, name = "grp_n"),
    by = "status"
  ) %>%
  mutate(hom_ref = grp_n - hom_alt - het) %>%
  select(variant, status, hom_ref, het, hom_alt) %>%
  pivot_wider(
    names_from  = status,
    values_from = c(hom_ref, het, hom_alt),
    names_sep   = "_",
    values_fill = 0
  ) %>%
  rename(
    d0 = hom_ref_case, d1 = het_case, d2 = hom_alt_case,
    h0 = hom_ref_control, h1 = het_control, h2 = hom_alt_control
  )

# Realizar test exacto alélico
results_allelic <- counts2 %>%
  rowwise() %>%
  mutate(p_allelic = allelic.exact.test(d0, d1, d2, h0, h1, h2)) %>%
  ungroup() %>%
  arrange(p_allelic)

# Seleccionar variantes significativas
sig_variants <- results_allelic %>%
  filter(p_allelic < 0.05)

variant_list <- sig_variants$variant

# Construir matriz binaria para ML (presencia/ausencia)
df_groups <- df_s %>% select(sample, patologia)

data_ml <- df_join %>%
  filter(sample %in% df_groups$sample, variant %in% variant_list) %>%
  mutate(pres = 1L) %>%
  distinct(sample, variant, pres) %>%
  complete(
    sample  = df_groups$sample,
    variant = variant_list,
    fill    = list(pres = 0L)
  ) %>%
  pivot_wider(names_from = variant, values_from = pres) %>%
  left_join(df_groups, by = "sample") %>%
  column_to_rownames("sample") %>%
  select(patologia, everything())

# Construir matriz de dosaje para RFE (0/1/2)
variant_list_all <- unique(df_variants$variant)
cases_info <- df_s %>% transmute(sample = as.character(sample), patologia)

data_rfe <- df_join %>%
  filter(status %in% c("case", "control")) %>%
  mutate(
    sample   = as.character(sample),
    variant  = as.character(variant),
    genotype = if_else(homo, 2L, 1L)
  ) %>%
  distinct(sample, variant, genotype) %>%
  complete(
    sample  = cases_info$sample,
    variant = variant_list_all,
    fill    = list(genotype = 0L)
  ) %>%
  pivot_wider(names_from = variant, values_from = genotype) %>%
  left_join(cases_info, by = "sample") %>%
  column_to_rownames("sample") %>%
  select(patologia, everything())

# Guardar datos para ML y RFE
saveRDS(data_ml,  "data_ml.rds")
saveRDS(data_rfe, "data_rfe.rds")

# Filtrar información de efectos y unir p-valores alélicos
effects_filtered <- df_effects %>%
  filter(id %in% variant_list) %>%
  arrange(desc(BIOTYPE == "miRNA")) %>%
  distinct(id, .keep_all = TRUE)

final_df <- effects_filtered %>%
  left_join(
    results_allelic %>% select(variant, p_allelic),
    by = c("id" = "variant")
  ) %>%
  select(
    id,
    SYMBOL,
    Consequence,
    VARIANT_CLASS,
    miRNA,
    pvalue = p_allelic
  ) %>%
  arrange(pvalue)

# Exportar resultados a archivo de texto
write.table(
  final_df,
  file      = "variantes_data.txt",
  sep       = "\t",
  row.names = FALSE,
  quote     = FALSE
)
