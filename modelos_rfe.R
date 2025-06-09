library(caret)
library(dplyr)
library(glmnet)
library(kernlab)
library(ranger)
library(doParallel)
library(pROC)
library(tidyr)
library(ggplot2)
library(stringr)
library(scales)

# Carga de datos y definición de la variable objetivo
data_rfe <- readRDS("data_rfe.rds")
data_rfe$patologia <- factor(
  data_rfe$patologia,
  levels = c("Hipertrofia", "Canalopatias")
)

# Función para balancear las clases mediante bootstrap 50/50
balancear_muestra <- function(df) {
  df_h <- df %>% filter(patologia == "Hipertrofia")
  df_c <- df %>% filter(patologia == "Canalopatias")
  n_c  <- nrow(df_c)
  samp_h <- df_h %>% slice_sample(n = n_c, replace = TRUE)
  bind_rows(samp_h, df_c) %>% slice_sample(prop = 1)
}

# Control de entrenamiento con validación cruzada repetida y down-sampling
ctrl <- trainControl(
  method          = "repeatedcv",
  number          = 10,
  repeats         = 3,
  summaryFunction = defaultSummary,
  classProbs      = TRUE,
  sampling        = "down",
  allowParallel   = TRUE
)

# Definición de las rejillas de hiperparámetros
grid_glmnet <- expand.grid(
  alpha  = seq(0, 1, length.out = 11),
  lambda = 10^seq(-5, 3, length.out = 50)
)

grid_svm <- expand.grid(
  C     = c(0.1, 1, 10),
  sigma = c(0.03619703, 0.1, 1)
)

p <- ncol(data_rfe) - 1
raiz_variables <- round(sqrt(p))
grid_rf <- expand.grid(
  mtry          = pmax(1, c(
    raiz_variables - 10,
    raiz_variables - 5,
    raiz_variables - 3,
    raiz_variables,
    raiz_variables + 3,
    raiz_variables + 6,
    raiz_variables + 10,
    raiz_variables + 20
  )),
  min.node.size = c(1, 2, 3, 4, 5),
  splitrule     = c("gini", "extratrees", "hellinger")
)

# Configuración de semillas, particiones y paralelismo
seeds  <- c(123, 456, 789)
splits <- c(0.70, 0.75, 0.80)
cores  <- parallel::detectCores(logical = FALSE) - 1
registerDoParallel(cores = cores)

# Carga de las variantes seleccionadas por RFE
rfe_vars <- readRDS("rfe_vars.rds")
rfe_vars <- lapply(rfe_vars, function(v) sub("^X", "", v))
rfe_methods <- names(rfe_vars)

# Entrenamiento y evaluación para cada método RFE
all_results <- list()

for (method in rfe_methods) {
  message("Método RFE: ", method)
  vars_sel <- rfe_vars[[method]]
  
  # Dataset con la patología y las variables seleccionadas
  df_sel <- data_rfe %>%
    select(patologia, all_of(vars_sel))
  
  for (seed in seeds) {
    message("  Semilla: ", seed)
    set.seed(seed)
    
    for (split in splits) {
      pct <- split * 100
      message("    Split: ", pct, "%")
      
      idx        <- createDataPartition(df_sel$patologia, p = split, list = FALSE)
      train_data <- df_sel[idx, ]
      test_data  <- df_sel[-idx, ]
      train_bal  <- balancear_muestra(train_data)
      
      specs <- list(
        LR  = list(method = "glmnet",   grid = grid_glmnet, pre = c("center", "scale")),
        SVM = list(method = "svmRadial", grid = grid_svm,    pre = c("center", "scale")),
        RF  = list(method = "ranger",    grid = grid_rf,     pre = NULL)
      )
      
      for (mn in names(specs)) {
        message("      Entrenando ", mn)
        t_start <- Sys.time()
        
        mod <- train(
          patologia ~ .,
          data       = train_bal,
          method     = specs[[mn]]$method,
          tuneGrid   = specs[[mn]]$grid,
          trControl  = ctrl,
          metric     = "Accuracy",
          preProcess = specs[[mn]]$pre,
          num.trees  = if (mn == "RF") 500 else NULL
        )
        
        train_time <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
        
        preds <- predict(mod, test_data)
        probs <- predict(mod, test_data, type = "prob")[, "Hipertrofia"]
        cm    <- confusionMatrix(preds, test_data$patologia, positive = "Hipertrofia")
        roc   <- roc(
          response  = test_data$patologia,
          predictor = probs,
          levels    = c("Canalopatias", "Hipertrofia"),
          direction = "<"
        )
        
        key <- paste(method, seed, pct, mn, sep = "_")
        all_results[[key]] <- list(
          model      = mod,
          train_time = train_time,
          metrics    = list(
            Accuracy    = as.numeric(cm$overall["Accuracy"]),
            Sensitivity = as.numeric(cm$byClass["Sensitivity"]),
            Specificity = as.numeric(cm$byClass["Specificity"]),
            AUC         = as.numeric(auc(roc))
          )
        )
        
        message(sprintf(
          "      %s: tiempo=%.1f s | Acc=%.3f | AUC=%.3f",
          mn, train_time, cm$overall["Accuracy"], auc(roc)
        ))
      }
    }
  }
}

# Guardar todos los resultados
saveRDS(all_results, "models_rfe.rds")
message("Resultados guardados en 'models_rfe.rds'")

# Preparación de los datos para el heatmap
res_list   <- readRDS("models_rfe.rds")
results_df <- bind_rows(
  lapply(names(res_list), function(key) {
    parts      <- strsplit(key, "_", fixed = TRUE)[[1]]
    rfe_method <- parts[1]
    seed       <- as.integer(parts[2])
    split_pct  <- as.numeric(parts[3])
    model_name <- parts[4]
    mets       <- res_list[[key]]$metrics
    tibble(
      RFE_Method  = rfe_method,
      Seed        = seed,
      Split       = split_pct,
      Model       = model_name,
      Accuracy    = mets$Accuracy,
      Sensitivity = mets$Sensitivity,
      Specificity = mets$Specificity
    )
  })
) %>%
  pivot_longer(
    cols      = c(Accuracy, Sensitivity, Specificity),
    names_to  = "Metric",
    values_to = "Value"
  ) %>%
  mutate(
    text_col = ifelse(Value > 0.5, "white", "black"),
    label    = sprintf("%.2f", Value)
  )

heat_df <- results_df %>%
  group_by(RFE_Method, Split, Model, Metric) %>%
  summarise(Value = mean(Value), .groups = "drop") %>%
  mutate(label = sprintf("%.2f", Value))

# Generación del heatmap de métricas
ggplot(heat_df, aes(x = Metric, y = factor(Split), fill = Value)) +
  geom_tile(color = "white") +
  geom_text(aes(label = label), color = "black", size = 4) +
  facet_grid(Model ~ RFE_Method, scales = "free") +
  scale_fill_gradient2(
    low      = "white",
    mid      = "white",
    high     = "red",
    midpoint = 0.5,
    limits   = c(0, 1),
    name     = "Media"
  ) +
  labs(
    x     = NULL,
    y     = "Split de entrenamiento (%)",
    title = "Métricas por Modelo, Split y Algoritmo RFE (media sobre 3 semillas)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1),
    panel.grid       = element_blank(),
    strip.background = element_rect(fill = "grey90", color = NA),
    strip.text.x     = element_text(face = "bold"),
    strip.text.y     = element_text(face = "bold"),
    plot.title       = element_text(face = "bold")
  )
