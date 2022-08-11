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

##Se remueven las presencias duplicadas
unic_list <- list()
for(i in 1:length(lista_especies)) {
    dupl <- duplicated(lista_especies[[i]][c("LONGITUD", "LATITUD")])
    pres_unic <- lista_especies[[i]][!dupl,]
    unic_list[[i]] <- pres_unic
}

##Datos erroneos (que tienen registro en el mar) son eliminados
for(i in 1:length(unic_list)) {
    unico <- unic_list[[i]][which(unic_list[[i]]$LATITUD >-110 &
                                      unic_list[[i]]$LONGITUD < -40),]
    condiciones <- extract(capasMXPresente, unico)
    regis_mal <- is.na(condiciones[,1])
    unic_list[[i]] <- unico[!regis_mal,]
}

##Se agrega sólo un dato por presencia (no es necesario tener más para generar
##los modelos)
pres_final <- list()
for(i in 1:length(unic_list)) {
    celdas <- cellFromXY(capasMXPresente[[1]], unic_list[[i]])
    dupl <- duplicated(celdas)
    pres_final[[i]] <- unic_list[[i]][!dupl,]
}

##Se generan datos "espaciales" a las ocurrencias, y se agrega el sistema de
##referencia de coordenadas
geo_pres <- list()
coord_ref <- CRS("+init=epsg:4326") 
for(i in 1:length(pres_final)) {
    unic <- pres_final[[i]]
    name <- as.data.frame(rep(especies_corr[[i]], 
                              nrow(unic)))
    spat_point <- SpatialPointsDataFrame(coords = unic, 
                                         data = name, 
                                         proj4string = coord_ref)
    geo_pres[[i]] <- spat_point
}

##Por especie, se guarda cerca del 30% de los datos para prueba del modelo
set.seed(1)
sp_train <- list()
sp_test <- list()
for (i in 1:length(geo_pres)) {
    df <- geo_pres[[i]]
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
##para los datos crudos es logistico. El multiplicador beta se ajusta a 0.3. Se
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

##Utilizando la evaluación de los modelos de entrenamiento, se calcula el umbral 
##de probabilidad para la clasificación binaria de presencia y ausencia para las
##proyecciones
umbral_train <- numeric()
umbral_test <- numeric()
umbral <- numeric()
for(i in 1:length(eval_train_list)) {
    umbraltrain_max <- eval_train_list[[i]]@t[which.max(
        eval_train_list[[i]]@TPR + eval_train_list[[i]]@TNR)]
    umbraltest_max <- eval_test_list[[i]]@t[which.max(
        eval_test_list[[i]]@TPR + eval_test_list[[i]]@TNR)]
    umbral_train[i] <- umbraltrain_max
    umbral_test[i] <- umbraltest_max
    umbral[i] <- (umbraltrain_max + umbraltest_max)/2
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
}
tabla1 <- data.frame("Especie" = especies_corr, 
                     "AUC entrenamiento" = AUC_train, 
                     "AUC prueba" = AUC_test)

##Se genera la tabla 1
tab_df(tabla1, 
       title = "Tabla 1. Lista de especies utilizadas en la generacion de modelos de distribucion de nicho, presentando los valores de AUC (area bajo la curva) para la evaluacion del rendimiento de los modelos. Se comparan los AUC para los sets de prueba y de entrenamiento del modelo.",
       alternate.rows = T) 

##El nombre de cada variable ambiental utilizada en el modelo de predicción
bio_var <- c("1 Temperatura media anual", 
             "2 Rango diurno medio", 
             "3 Isotermalidad",
             "4 Estacionalidad de temperatura (desviación estándar)",
             "5 Temperatura máxima del mes más cálido",
             "6 Temperatrua mínima del mes más frío",
             "7 Rango anual de la temperatura",
             "8 Temperatura media del cuarto más húmedo",
             "9 Temperatura media del cuarto más seco",
             "10 Temperatura media del cuarto más cálido",
             "11 Temperatura media del cuarto más frio",
             "12 Precipitación anual",
             "13 Precipitación del mes más húmedo",
             "14 Precipitación del mes más seco",
             "15 Estacionalidad de temperatura (coeficiente de variación)",
             "16 Precipitación del cuarto más húmedo",
             "17Precipitación del cuarto más seco",
             "18 Precipitación del cuarto más cálido",
             "19 Precipitación del cuarto más frío")

##Se obtiene la importancia para cada variable entre todos los modelos
capas_contr <- data.frame("Capa" = c(1:19))
for(i in 1:length(modelo_sp)) {
    contr <- modelo_sp[[i]]@results
    rows <- grep(row.names(contr), pattern = ".contribution")
    capas_contr[especies_corr[i]] <- contr[c(rows),]
}
capas_contr <- capas_contr[, -1]

##Los promedios para la contribución de cada variable es calculada y tabulada
means_var <- numeric()
for(i in 1:nrow(capas_contr)) {
    mean_c <- rowMeans(capas_contr[i,])
    means_var[i] <- mean_c
}
mean_contr <- data.frame("Capa" = c("Bio 1", "Bio 2", "Bio 3", "Bio 4", "Bio 5", 
                                    "Bio 6", "Bio 7", "Bio 8", "Bio 9", "Bio 10",
                                    "Bio 11", "Bio 12", "Bio 13", "Bio 14", 
                                    "Bio 15", "Bio 16", "Bio 17", "Bio 18",
                                    "Bio 19"), 
                         "Promedio por modelo" = means_var)

##Se genera la gráfica 1
graph1 <- ggplot(data = mean_contr, 
            aes(y = reorder(Capa, Promedio.por.modelo), 
                x = Promedio.por.modelo))
graph1 + 
    geom_col() + 
    theme_bw() + 
    coord_cartesian(xlim = c(1.1, 26)) +
    labs(x = "Media de contribución para modelos", y = "Capa ambiental")

##Se obtienen los valores más altos y más bajos de contribución para cada capa
##ambiental y se tabulan
max_val <- numeric()
min_val <- numeric()
for(i in 1:nrow(capas_contr)) {
    x_max <- max(capas_contr[i,])
    x_min <- min(capas_contr[i,])
    max_val[i] <- x_max
    min_val[i] <- x_min
}
tabla2 <- data.frame("Capa" = c("Bio 1", "Bio 2", "Bio 3", "Bio 4", "Bio 5", 
                                 "Bio 6", "Bio 7", "Bio 8", "Bio 9", "Bio 10",
                                 "Bio 11", "Bio 12", "Bio 13", "Bio 14", 
                                 "Bio 15", "Bio 16", "Bio 17", "Bio 18",
                                 "Bio 19"), 
                      "Minimo" = min_val,
                      "Maximo" = max_val) 
tab_df(tabla2, 
       alternate.rows = T) 
