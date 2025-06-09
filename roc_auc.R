library(dplyr)
library(ggplot2)
library(pROC)
library(stringr)
library(scales)
library(caret)

# Leer resultados y datos
res     <- readRDS("models_results.rds")
data_ml <- readRDS("data_ml.rds")
data_ml$patologia <- factor(
  data_ml$patologia,
  levels = c("Canalopatias", "Hipertrofia")
)

seeds  <- c(123, 456, 789)
splits <- c(0.70, 0.75, 0.80)
models <- c("LR", "SVM", "RF")

# Inicializar lista para guardar predicciones
all_preds <- lapply(models, function(m) {
  sapply(c("70", "75", "80"), function(pct) {
    list(truth = character(), probs = numeric())
  }, simplify = FALSE)
})
names(all_preds) <- models

# Reproducir particiones, predecir y acumular resultados
for (seed in seeds) {
  set.seed(seed)
  for (split in splits) {
    split_pct <- as.integer(split * 100)
    idx       <- createDataPartition(data_ml$patologia, p = split, list = FALSE)
    test_data <- data_ml[-idx, ]
    
    for (mn in models) {
      key  <- paste(seed, split, mn, sep = "_")
      mod  <- res$models[[key]]
      
      probs <- tryCatch({
        predict(mod, test_data, type = "prob")[, "Hipertrofia"]
      }, error = function(e) {
        warning("Predicción fallida para ", key, ": ", e$message)
        rep(NA_real_, nrow(test_data))
      })
      if (all(is.na(probs))) next
      
      truth <- as.character(test_data$patologia)
      
      all_preds[[mn]][[as.character(split_pct)]]$probs <- 
        c(all_preds[[mn]][[as.character(split_pct)]]$probs, probs)
      all_preds[[mn]][[as.character(split_pct)]]$truth <- 
        c(all_preds[[mn]][[as.character(split_pct)]]$truth, truth)
    }
  }
}

# Asegurar niveles consistentes en las etiquetas verdaderas
all_preds <- lapply(all_preds, function(splits_list) {
  lapply(splits_list, function(entry) {
    entry$truth <- factor(entry$truth, levels = c("Canalopatias", "Hipertrofia"))
    entry
  })
})

# Calcular curvas ROC y AUC para cada modelo y split
roc_df_list <- list()
auc_table    <- data.frame(
  Model = character(),
  Split = integer(),
  AUC   = numeric(),
  stringsAsFactors = FALSE
)

for (mn in models) {
  for (split_pct in c("70", "75", "80")) {
    entry  <- all_preds[[mn]][[split_pct]]
    truths <- entry$truth
    probs  <- entry$probs
    
    if (length(unique(truths)) < 2) {
      message("Omitiendo ROC para ", mn, " split=", split_pct, " (una sola clase)")
      next
    }
    
    roc_obj <- roc(
      response  = truths,
      predictor = probs,
      levels    = c("Canalopatias", "Hipertrofia"),
      direction = "<"
    )
    auc_val <- as.numeric(roc_obj$auc)
    
    auc_table <- rbind(
      auc_table,
      data.frame(
        Model = mn,
        Split = as.integer(split_pct),
        AUC   = auc_val,
        stringsAsFactors = FALSE
      )
    )
    
    df_temp <- data.frame(
      Model = mn,
      Split = as.integer(split_pct),
      fpr   = 1 - roc_obj$specificities,
      tpr   = roc_obj$sensitivities
    )
    roc_df_list[[paste(mn, split_pct, sep = "_")]] <- df_temp
  }
}

# Unir todos los data.frames de ROC
roc_df <- bind_rows(roc_df_list)

# Añadir columna combinada para etiquetas
roc_df   <- roc_df   %>% mutate(Model_Split = paste0(Model, "_", Split))
auc_table <- auc_table %>% mutate(Model_Split = paste0(Model, "_", Split))

# Crear mapa de etiquetas para la leyenda
# (se asume que 'Label' ya está definido en auc_table)
label_map <- setNames(auc_table$Label, auc_table$Model_Split)

# Dibujar curvas ROC
ggplot(roc_df, aes(x = fpr, y = tpr, color = Model_Split, group = Model_Split)) +
  geom_line(size = 1) +
  geom_abline(linetype = "dashed", color = "grey50") +
  scale_color_manual(
    name   = "ROC–AUC",
    values = c(
      "LR_70"  = "#D73027", "LR_75"  = "#F46D43", "LR_80"  = "#FDAE61",
      "SVM_70" = "#4575B4", "SVM_75" = "#74ADD1", "SVM_80" = "#ABD9E9",
      "RF_70"  = "#006400", "RF_75"  = "#31A354", "RF_80"  = "#A6D96A"
    ),
    labels = label_map,
    guide  = guide_legend(
      ncol   = 2,
      title.position = "top",
      byrow  = FALSE
    )
  ) +
  labs(
    title = "Curvas ROC por Modelo y Split",
    x     = "False Positive Rate",
    y     = "True Positive Rate"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position      = c(0.95, 0.05),
    legend.justification = c(1, 0),
    plot.title           = element_text(face = "bold", hjust = 0.5)
  )
