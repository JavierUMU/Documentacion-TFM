library(caret)
library(dplyr)
library(doParallel)
library(randomForest)  
library(ipred)         
library(klaR)          
library(pROC)

# Leer y preparar los datos
data_rfe <- readRDS("data_rfe.rds") %>%
  mutate(patologia = factor(patologia, levels = c("Hipertrofia", "Canalopatias")))

# Asegurar nombres válidos y convertir a data.frame
names(data_rfe) <- make.names(names(data_rfe), unique = TRUE)
data_rfe <- as.data.frame(data_rfe)

# Dividir en entrenamiento (80%) y prueba (20%)
set.seed(123)
idx       <- createDataPartition(data_rfe$patologia, p = 0.80, list = FALSE)
train_all <- data_rfe[idx, ]
test_all  <- data_rfe[-idx, ]

# Separar clases para balanceo interno
canalopatias <- filter(train_all, patologia == "Canalopatias")
hipertrofia  <- filter(train_all, patologia == "Hipertrofia")
n_control    <- nrow(canalopatias)

# Definir métodos RFE y sus funciones
methods    <- c("rf", "treebag", "nb")
funcs_list <- list(
  rf      = rfFuncs,
  treebag = treebagFuncs,
  nb      = nbFuncs
)

# Contenedores para guardar resultados
final_vars    <- list()
final_metrics <- tibble(
  Method      = character(),
  Accuracy    = numeric(),
  Sensitivity = numeric(),
  Specificity = numeric(),
  AUC         = numeric()
)

# Aplicar RFE balanceado y evaluar cada método
for (m in methods) {
  message("Método RFE: ", m)
  
  # Configuración de RFE con validación cruzada
  ctrl <- rfeControl(
    functions     = funcs_list[[m]],
    method        = "cv",
    number        = 5,
    allowParallel = TRUE
  )
  
  # Repetir RFE varias veces para identificar variables estables
  sel_runs <- vector("list", 5)
  for (i in seq_len(5)) {
    it0    <- Sys.time()
    samp_h <- slice_sample(hipertrofia, n = n_control)
    data_bal <- bind_rows(samp_h, canalopatias)
    x_rfe <- as.data.frame(select(data_bal, -patologia))
    y_rfe <- data_bal$patologia
    sizes <- seq(5, ncol(x_rfe), by = 5)
    
    set.seed(2000 + i)
    rfe_res <- tryCatch(
      rfe(x = x_rfe, y = y_rfe, sizes = sizes, rfeControl = ctrl),
      error = function(e) {
        message("  Error RFE iter ", i, ": ", e$message)
        NULL
      }
    )
    
    sel_runs[[i]] <- if (!is.null(rfe_res)) rfe_res$optVariables else character(0)
    iti <- as.numeric(difftime(Sys.time(), it0, units = "secs"))
    message("  Iteración ", i, " en ", round(iti, 1), " s")
  }
  
  # Seleccionar las variables que aparecen al menos en 3 iteraciones
  freqs <- table(unlist(sel_runs))
  final_vars[[m]] <- names(freqs)[freqs >= 3]
  message("Variables finales (", m, "): ", length(final_vars[[m]]))
  
  # Entrenar modelo final con las variables seleccionadas
  vars_sel  <- final_vars[[m]]
  train_fit <- select(train_all, patologia, all_of(vars_sel))
  
  set.seed(3210)
  model <- train(
    patologia ~ .,
    data       = train_fit,
    method     = switch(m, rf = "ranger", treebag = "treebag", nb = "nb"),
    trControl  = trainControl(classProbs = TRUE),
    tuneLength = 3
  )
  
  # Evaluar el modelo con el conjunto de prueba
  test_fit <- select(test_all, patologia, all_of(vars_sel))
  preds <- predict(model, test_fit)
  probs <- predict(model, test_fit, type = "prob")[, "Hipertrofia"]
  
  cm  <- confusionMatrix(preds, test_fit$patologia, positive = "Hipertrofia")
  roc <- roc(
    response  = test_fit$patologia,
    predictor = probs,
    levels    = c("Canalopatias", "Hipertrofia"),
    direction = "<"
  )
  
  # Guardar métricas finales
  final_metrics <- final_metrics %>%
    add_row(
      Method      = m,
      Accuracy    = as.numeric(cm$overall["Accuracy"]),
      Sensitivity = as.numeric(cm$byClass["Sensitivity"]),
      Specificity = as.numeric(cm$byClass["Specificity"]),
      AUC         = as.numeric(auc(roc))
    )
  
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  message("Modelo ", m, " completado en ", round(elapsed, 1), " min")
}

# Guardar variables seleccionadas
saveRDS(final_vars, "rfe_vars.rds")

# Finalizar cluster paralelo
stopCluster(cl)

# Mostrar resumen de número de variables por método
df_summary <- tibble(
  Method    = names(final_vars),
  Variables = lengths(final_vars)
)

print(df_summary)
