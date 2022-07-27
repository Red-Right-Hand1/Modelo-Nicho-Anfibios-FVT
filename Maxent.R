##LIBRERIAS
library(dismo)
library(raster)
library(plyr)
library(rJava)
library(sf)
library(rgeos)
library(rgdal)
library(sjPlot)
library(stargazer)

###LIMPIEZA DE DATOS

##Directorio de trabajo
setwd("C:/Users/Tabby/Documents/Modelos")

##Se crea una función para la extracción de los datos Raster para cada region
##y cada tiempo. Arroja una lista con los datos Raster
rasUpload <- function(region, tiempo){ ##region = "Eje" o "MX"; tiempo = "PRESENTE" o "FUTURO"
  wd = paste0(
    "C:/Users/Tabby/Documents/Modelos/Capas/", 
    tiempo, "/", region, "_30s/ASC") ##Se genera el directorio para buscar
  raster_name <- grep(
    list.files(path = wd), 
    pattern = ".asc", 
    value = TRUE) ##Se buscan los nombres de los archivos asc (que ya tienen sus proyecciones)
  capas <- list() ##Lista donde se van a almacenar los datos Raster
  for(i  in 1:length(raster_name)){ ##loop para buscar y convertir a Raster cada archivo
    Ras <- raster(paste0(wd, "/", raster_name[i]))    
    capas[[i]] <- Ras
  } 
  return(capas) ##La función devuelve la lista con los archivos Raster
}

##Se crea la lista para el Eje en el presente
rasterEjePresente <- rasUpload("Eje", "PRESENTE")

##Se crea la lista para el Eje en el futuro
rasterEjeFuturo <- rasUpload("Eje", "FUTURO")

##Se crea la lista para MX en el presente
rasterMXPresente <- rasUpload("MX", "PRESENTE")

##Posteriormente se apilan todas las capas ambientales por región y tiempo
capasEjePresente <- stack(rasterEjePresente)
capasEjeFuturo <- stack(rasterEjeFuturo)
capasMXPresente <- stack(rasterMXPresente)

##Se delimitan los puntos de fondo
set.seed(1)
bp <- sampleRandom(x = capasMXPresente,
                   size = 100000,
                   na.rm = T,
                   sp = T)

##Se carga el archivo csv de ocurrencias
ocurrencias <- read.csv(
  "C:/Users/Tabby/Documents/Modelos/Ocurrencias/143_MX.xlsx - Hoja1.csv")

##Un data frame para cada especie es creada
especies <- count(ocurrencias[,1])[c(47:91),1]
especies_corr <- grep(
  especies, pattern = 
    "Eleutherodactylus dilatus|Eleutherodactylus saxatilis|Isthmura naucampatepetl|Eleutherodactylus grandis|Eleutherodactylus maurus",
  invert = TRUE, value = TRUE)
lista_especies = list()
for(i in 1:length(especies_corr)){
  numero_fila <- which(ocurrencias$ESPECIE == especies_corr[i])
  df <- ocurrencias[c(numero_fila[1]: numero_fila[length(numero_fila)]), c(2,3)]
  lista_especies[[i]] <- df
}
names(lista_especies) <- especies_corr

##Por especie, se guarda cerca del 30% de los datos para prueba del modelo
set.seed(1)
sp_train <- list()
sp_test <- list()
for (i in 1:length(lista_especies)) {
  df <- lista_especies[[i]]
  selected <- sample(1:nrow(df), nrow(df)*0.3)
  df_train <- df[-selected,]
  df_test <- df[selected,]
  sp_train[[i]] <- df_train
  sp_test[[i]] <- df_test
}
names(sp_train) <- especies_corr
names(sp_test) <- especies_corr

##Se extraen las condiciones ambientales para los datos de ocurrencia de 
##entrenamiento y de prueba, para cada especie
p_list <- list()
p_test_list <- list()
for(i in 1:length(lista_especies)) {
  p <- as.data.frame(extract(capasMXPresente, sp_train[[i]]))
  p_test <- as.data.frame(extract(capasMXPresente, sp_test[[i]]))
  p_list[[i]] <- p
  p_test_list[[i]] <- p_test
}
names(p_list) <- especies_corr
names(p_test_list) <- especies_corr

##De igual manera, se extraen las condiciones ambientales del fondo
a <- as.data.frame(extract(capasMXPresente, bp))

##Los datos para MaxEnt son transformados en formato tabular (binario, usando
##1 = presencia y 0 = pseudo-ausencia), asignando 1 a los datos ambientales 
##obtenidos de los datos de entrenamiento, mientras que 0 son asignados a 
##los datos ambientales obtenidos del fondo. De igual manera, se crea un data 
##frame utilizando los datos tabulados de presencia con los datos de condiciones 
##ambientales de fondo
pa_list <- list()
pder_list <- list()
for(i in 1:length(lista_especies)) {
  pres_aus <- c(rep(1, nrow(p_list[[i]])), rep(0, nrow(a)))
  pa_list[[i]] <- pres_aus 
  pder_list[[i]] <- as.data.frame(rbind(p_list[[i]], a))
} 
names(pa_list) <- especies_corr ##Datos de ausencia/presencia
names(pder_list) <- especies_corr ##Valores ambientales correspondientes a los datos de ausencia/presencia

###MODELAJE Y PREDICCIÓN DE DISTRIBUCIÓN

##La función para modificar los parámetros de MaxEnt es cargada (esta función
##se puede obtener por medio del repositorio de la clase de "R and Maxent" de
##NIMBioS: https://github.com/shandongfx/nimbios_enm)
setwd("C:/Users/Tabby/Documents/Modelos")
source("Appendix2_prepPara.R")

##Se genera el modelo calibrado a las capas ambientales de México al presente, 
##para cada especie. Los parámetros para el modelaje toman en cuenta restricciones
##como las Lineales y Cuadráticas. Se aplica jackknife y el formato de salida 
##para los datos crudos es logistico. El multiplicador beta se ajusta a 0.4. Se
##aplica clamp, crossvalidate y extrapolación
wd <- "C:/Users/Tabby/Documents/Modelos/Modelo/"
model_wd <- paste0(wd, as.character(length(list.files(wd)) + 1), "/")
modelo_sp <- list()
for(i in 1:length(sp_train)){
  modelo <- maxent(
    x = pder_list[[i]],
    p = pa_list[[i]],
    path = paste0(model_wd,
                  especies_corr[i]),
    removeDuplicates = TRUE,
    args = prepPara(userfeatures = "LQ", 
                    jackknife = TRUE,
                    outputfiletype = "asc",
                    outputformat = "logistic",
                    betamultiplier = 0.3, 
                    randomseed = TRUE, 
                    doclamp = TRUE, 
                    replicatetype = "crossvalidate",
                    extrapolate = TRUE)
  )
  modelo_sp[[i]] <- modelo
}
names(modelo_sp) <- especies_corr

##Se genera la proyección de distribución por especie hacia las capas 
##abioticas del Eje en el presente y para el futuro
proyec_presente <- list()
proyec_futuro <- list()
for (i in 1:length(modelo_sp)) {
  proyeccion <- predict(modelo_sp[[i]], capasEjePresente)
  proyec_presente[[i]] <- proyeccion
}
for (i in 1:length(modelo_sp)) {
  proyeccion <- predict(modelo_sp[[i]], capasEjeFuturo)
  proyec_futuro[[i]] <- proyeccion
}

##Los modelos son evaluados, respecto al set de entrenamiento y de prueba
eval_train_list <- list()
eval_test_list <- list()
for(i in 1:length(modelo_sp)) {
  eval_train <- evaluate(p = p_list[[i]], a = a, model = modelo_sp[[i]])
  eval_test <- evaluate(p = p_test_list[[i]], a = a, model = modelo_sp[[i]])
  eval_train_list[[i]] <- eval_train
  eval_test_list[[i]] <- eval_test
}

##Los mapas de predicción son guardados en las carpetas correspondientes
wd_presente <- "C:/Users/Tabby/Documents/Modelos/Proyecciones/Presente/"
wd_futuro <- "C:/Users/Tabby/Documents/Modelos/Proyecciones/Fututo/"
proy_presente_wd <- paste0(wd_presente,
                           as.character(length(list.files(wd_presente)) + 1), "/")
proy_futuro_wd <- paste0(wd_futuro,
                         as.character(length(list.files(wd_futuro)) + 1), "/")
dir.create(proy_presente_wd)
dir.create(proy_futuro_wd)
for(i in 1:length(proyec_presente)){
  dir.create(paste0(proy_presente_wd, especies_corr[i]))
  writeRaster(proyec_presente[[i]],
              filename = paste0(proy_presente_wd, especies_corr[i], 
                                "/response_", 
                                especies_corr[i], ".asc"),
              format = "ascii",
              bylayer = TRUE,
              overwrite = T)
  png(file = paste0(proy_presente_wd, especies_corr[i], "/response_", 
                    especies_corr[i], ".png"),
      width = 1500, height = 1000)
  plot(proyec_presente[[i]])
  dir.create(paste0(proy_futuro_wd, especies_corr[i]))
  writeRaster(proyec_futuro[[i]],
              filename = paste0(proy_futuro_wd, especies_corr[i], 
                                "/response_", 
                                especies_corr[i], ".asc"),
              format = "ascii",
              bylayer = TRUE,
              overwrite = T)
  png(file = paste0(proy_futuro_wd, especies_corr[i], "/response_", 
                    especies_corr[i], ".png"),
      width = 1500, height = 1000)
  plot(proyec_futuro[[i]])
  dev.off()
}

###GENERACIÓN DE TABLAS Y GRÁFICAS

##Obtenemos los valores de AUC para cada evaluación del modelo, para cada
##especie. Agregamos esto a un data frame con el nombre de la especie 
##correspondiente
AUC_train <- integer()
AUC_test <- integer()
for(i in 1:length(eval_train_list)) {
  AUCtr <- eval_train_list[[i]]@auc
  AUCte <- eval_test_list[[i]]@auc
  AUC_train[i] <- AUCtr
  AUC_test[i] <- AUCte
  pa
}
tabla1 <- data.frame("Especie" = especies_corr, 
                     "AUC entrenamiento" = AUC_train, 
                     "AUC prueba" = AUC_test)

##Se genera la tabla 1
tab_df(tabla1, 
       title = "Tabla 1. Lista de especies utilizadas en la generación de modelos de distribución de nicho, presentando los valores de AUC (área bajo la curva) para la evaluación del rendimiento de los modelos. Se comparan los AUC para los sets de prueba y de entrenamiento del modelo.",
       alternate.rows = T,
       file = "tabla1.html") 

##El nombre de cada variable ambiental utilizada en el modelo de predicción
bio_var <- c("Temperatura media anual", 
             "Rango diurno medio", 
             "Isotermalidad",
             "Estacionalidad de temperatura (desviación estándar)",
             "Temperatura máxima del mes más cálido",
             "Temperatrua mínima del mes más frío",
             "Rango anual de la temperatura",
             "Temperatura media del cuarto más húmedo",
             "Temperatura media del cuarto más seco",
             "Temperatura media del cuarto más cálido",
             "Temperatura media del cuarto más frio",
             "Precipitación anual",
             "Precipitación del mes más húmedo",
             "Precipitación del mes más seco",
             "Estacionalidad de temperatura (coeficiente de variación)",
             "Precipitación del cuarto más húmedo",
             "Precipitación del cuarto más seco",
             "Precipitación del cuarto más cálido",
             "Precipitación del cuarto más frío")

##Se obtiene la importancia para cada variable entre todos los modelos
modelo_sp[[1]]

proyec_presente[[1]]
