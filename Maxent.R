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
        "Eleutherodactylus dilatus|Eleutherodactylus saxatilis|Isthmura naucampatepetl|Eleutherodactylus grandis|Eleutherodactylus maurus", invert = TRUE, value = TRUE)
lista_especies = list()
for(i in 1:length(especies_corr)){
    numero_fila <- which(ocurrencias$ESPECIE == especies_corr[i])
    df <- ocurrencias[c(numero_fila[1]: numero_fila[length(numero_fila)]), c(2,3)]
    lista_especies[[i]] <- df
}
names(lista_especies) <- especies_corr

##Por especie, se guarda cerca del 30% de los datos para prueba del modelo




dfprueba <- lista_especies[[1]]
fold <- kfold(dfprueba, k = 3.4)
occ_train <- dfprueba[fold != 1,]
occ_test <- dfprueba[fold == 1,]

me <- maxent(
    x = capasMXPresente, 
    p = occ_train,
    nbp = 100000,
    removeDuplicates = TRUE)
plot(me)
response(me)

r <- predict(me, capasEjePresente)
plot(r)
points(dfprueba)

##Test
bp <- randomPoints(capasEjePresente, 10000)

el <- evaluate(me, p = occ_test, a = bp, x = capasEjePresente)
el




###MODELAJE Y PREDICCIÓN DE DISTRIBUCIÓN

##Se genera el modelo calibrado a las capas ambientales de México al presente, 
##para cada especie
modelo <- maxent(
    x = capasMXPresente, 
    p = lista_especies$`Dermophis oaxacae`, 
    nbp = 100000, 
    removeDuplicates = TRUE,
    path = "C:/Users/Tabby/Documents/Modelos/Proyecciones/Presente/R/01")




##Preguntar sobre la lista de capas