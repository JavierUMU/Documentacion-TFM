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

## Cargar datos y definir variable de interés como factor
data_ml <- readRDS("data_ml.rds")
data_ml$patologia <- factor(data_ml$patologia,
                            levels = c("Hipertrofia", "Canalopatias"))

## Función para balancear la muestra (50/50) mediante bootstrap del grupo mayoritario
balancear_muestra <- function(df) {
  df_h <- df %>% filter(patologia == "Hipertrofia")
  df_c <- df %>% filter(patologia == "Canalopatias")
  n_c  <- nrow(df_c)
  samp_h <- df_h %>% slice_sample(n = n_c, replace = TRUE)
  bind_rows(samp_h, df_c) %>% slice_sample(prop = 1)
}

## Configurar validación cruzada repetida con down-sampling y cálculo de AUC
ctrl <- trainControl(
  method          = "repeatedcv",
  number          = 10,
  repeats         = 3,
  summaryFunction = twoClassSummary,  # métricas ROC/AUC
  classProbs      = TRUE,             # obtener probabilidades
  sampling        = "down",
  allowParallel   = TRUE,
  verboseIter     = FALSE             # sin mostrar iteraciones
)

## Definir cuadrículas de hiperparámetros para cada modelo
grid_glmnet <- expand.grid(
  alpha  = seq(0, 1, length.out = 11),
  lambda = 10^seq(-5, 3, length.out = 50)
)
grid_svm <- expand.grid(
  C     = c(0.1, 1, 10),
  sigma = c(0.03619703, 0.1, 1)
)

p <- ncol(data_ml) - 1
raiz_variables <- round(sqrt(p))
grid_rf <- expand.grid(
  mtry          = pmax(1, raiz_variables + c(-10, -5, -3, 0, 3, 6, 10, 20)),
  min.node.size = c(1, 2, 3, 4, 5),
  splitrule     = c("gini", "extratrees", "hellinger")
)

## Definir semillas y proporciones de split para experimentos
seeds  <- c(123, 456, 789)
splits <- c(0.70, 0.75, 0.80)

## Preparar estructuras para almacenar modelos, métricas y tiempos
all_results <- list(
  models  = list(),
  metrics = tibble(
    Seed        = integer(),
    Split       = numeric(),
    Model       = character(),
    Accuracy    = numeric(),
    Sensitivity = numeric(),
    Specificity = numeric(),
    AUC         = numeric()
  ),
  times   = tibble(
    Seed      = integer(),
    Split     = numeric(),
    Model     = character(),
    TrainTime = numeric()
  )
)

## Bucle principal: entrenar y evaluar cada modelo para cada semilla y split
for (seed in seeds) {
  set.seed(seed)
  for (split in splits) {
    message("Seed=", seed, " Split=", split*100, "%")
    idx        <- createDataPartition(data_ml$patologia, p = split, list = FALSE)
    train_data <- data_ml[idx, ]
    test_data  <- data_ml[-idx, ]
    train_bal  <- balancear_muestra(train_data)
    
    specs <- list(
      LR  = list(method = "glmnet",   grid = grid_glmnet, pre = c("center", "scale")),
      SVM = list(method = "svmRadial", grid = grid_svm,   pre = c("center", "scale")),
      RF  = list(method = "ranger",    grid = grid_rf,    pre = NULL)
    )
    
    for (mn in names(specs)) {
      message(" -> Entrenando ", mn)
      t0 <- Sys.time()
      mod <- suppressWarnings(
        suppressMessages(
          train(
            patologia ~ .,
            data       = train_bal,
            method     = specs[[mn]]$method,
            tuneGrid   = specs[[mn]]$grid,
            trControl  = ctrl,
            metric     = "Accuracy",
            preProcess = specs[[mn]]$pre,
            num.trees  = if (mn == "RF") 500 else NULL
          )
        )
      )
      secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      message(" <- ", mn, " en ", round(secs, 1), " s")
      
      # Guardar modelo y tiempo de entrenamiento
      key <- paste(seed, split, mn, sep = "_")
      all_results$models[[key]] <- mod
      all_results$times <- add_row(
        all_results$times,
        Seed      = seed,
        Split     = split * 100,
        Model     = mn,
        TrainTime = secs
      )
      
      # Predecir en datos de test y calcular métricas
      pred   <- factor(predict(mod, test_data), levels = levels(test_data$patologia))
      probs  <- predict(mod, test_data, type = "prob")[, "Hipertrofia"]
      cm     <- confusionMatrix(pred, test_data$patologia, positive = "Hipertrofia")
      roc_obj <- roc(
        response  = test_data$patologia,
        predictor = probs,
        levels    = c("Canalopatias", "Hipertrofia")
      )
      auc_val <- as.numeric(auc(roc_obj))
      
      # Guardar métricas en la tabla
      all_results$metrics <- add_row(
        all_results$metrics,
        Seed        = seed,
        Split       = split * 100,
        Model       = mn,
        Accuracy    = cm$overall["Accuracy"],
        Sensitivity = cm$byClass["Sensitivity"],
        Specificity = cm$byClass["Specificity"],
        AUC         = auc_val
      )
    }
  }
}

## Guardar resultados en disco
saveRDS(all_results, "models_results.rds")
message("Resultados guardados en 'models_results.rds'")

## Generar heatmap de métricas promedio por modelo y split
# Cargar resultados y datos originales
res     <- readRDS("models_results.rds")
# Calcular medias por combinación de split y modelo
metrics_avg <- res$metrics %>%
  group_by(Split, Model) %>%
  summarise(
    Accuracy    = mean(Accuracy),
    Sensitivity = mean(Sensitivity),
    Specificity = mean(Specificity),
    .groups     = "drop"
  )

# Convertir a formato largo para ggplot
metrics_long <- metrics_avg %>%
  pivot_longer(
    cols      = Accuracy:Specificity,
    names_to  = "Metric",
    values_to = "MeanValue"
  )

# Dibujar heatmap con etiquetas de valor
ggplot(metrics_long, aes(x = Metric, y = factor(Split), fill = MeanValue)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", MeanValue)), size = 3) +
  facet_grid(. ~ Model) +
  scale_fill_gradient(low = "white", high = "red", name = "Media") +
  labs(
    x     = NULL,
    y     = "Split de entrenamiento (%)"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid       = element_blank(),
    axis.text.x      = element_text(angle = 45, hjust = 1),
    strip.background = element_rect(fill = "grey90", color = NA),
    strip.text       = element_text(face = "bold")
  )
