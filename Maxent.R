##Librerías
library(dismo)
library(raster)
library(plyr)
library(rJava)

###LIMPIEZA DE DATOS

##Directorio de trabajo
setwd("C:/Users/Tabby/Documents/Modelos")

##Se crea una función para la extracción de los datos Raster para cada region
##y cada tiempo. Arroja una lista con los datos Raster
rasUpload <- function(region, tiempo){ ##region = "Eje" o "MX"; tiempo = "PRESENTE" o "FUTURO"
    wd = paste0("C:/Users/Tabby/Documents/Modelos/Capas/", tiempo, "/", region, "_30s") ##Se genera el directorio para buscar
    exclusion <- grep(list.files(path = wd), pattern = ".aux.xml|.ovr", invert = TRUE, value = TRUE) ##Se excluyen los archivos con terminación ".aux.xml" y ".ovr"
    raster_name <- grep(exclusion, pattern = "bio_", value = TRUE) ##Se buscan los nombres de los folders con los archivos raster
    capas <- list() ##Lista donde se van a almacenar los datos Raster
    for(i  in 1:length(raster_name)){ ##loop para buscar y convertir a Raster cada folder
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

##Se carga el archivo csv de ocurrencias
ocurrencias <- read.csv("C:/Users/Tabby/Documents/Modelos/Ocurrencias/143_MX.xlsx - Hoja1.csv")

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
sp_train <- list()
sp_test <- list()
for (i in 1:length(lista_especies)) {
    df <- lista_especies[[i]]
    fold <- kfold(df, k = 3.4)
    df_train <- df[fold != 1,]
    df_test <- df[fold == 1,]
    sp_train[[i]] <- df_train
    sp_test[[i]] <- df_test
}
names(sp_train) <- especies_corr
names(sp_test) <- especies_corr

###MODELAJE Y PREDICCIÓN DE DISTRIBUCIÓN

##Se genera el modelo calibrado a las capas ambientales de México al presente, 
##para cada especie
modelo_sp <- list()
for(i in 1:length(sp_train)){
    modelo <- maxent(
        x = capasMXPresente,
        p = sp_train[[i]],
        nbp = 100000,
        removeDuplicates = TRUE,
        path = paste0("C:/Users/Tabby/Documents/Modelos/Proyecciones/Presente/R/01/", especies_corr[i])
    )
    modelo_sp[[i]] <- modelo
}
names(modelo_sp) <- especies_corr

##Se grafican las curvas de respuesta para cada especie, y se almacena 
##en su carpeta correspondiente
for(i in 1:length(modelo_sp)){
    png(file = paste0("C:/Users/Tabby/Documents/Modelos/Proyecciones/Presente/R/01/", especies_corr[i], "/response_", especies_corr[i], ".png"),
        width=1500, height=1000)
    response(modelo_sp[[i]])
    dev.off()
}

##Se genera la proyección de distribución por especie hacia las capas 
##abioticas del Eje al presente
proyec_presente <- list()
for (i in 1:length(modelo_sp)) {
    proyeccion <- predict(modelo_sp[[i]], capasEjePresente)
    proyec_presente[[i]] <- proyeccion
}

##Los mapas de predicción son guardados en las carpetas correspondientes
for(i in 1:length(proyec_presente)){
    png(file = paste0("C:/Users/Tabby/Documents/Modelos/Proyecciones/Presente/R/01/", especies_corr[i], "/response_", especies_corr[i], ".png"),
        width=1500, height=1000)
    plot(proyec_presente[[i]])
    points(proyec_presente[[i]])
    dev.off()
}
