
###### 1. Instalación de paquetes necesarios para la realización de este ejercicio#####

# Geo Query
# http://www.bioconductor.org/packages/release/bioc/html/GEOquery.html
# Paquete que permite descargar del repositorio GEO (Gene Expression Omnibus) conjuntos de microarrays
if (!("GEOquery" %in% installed.packages())) { 
  source("http://www.bioconductor.org/biocLite.R"); 
  biocLite("GEOquery");
}


# AFFY
# http://www.bioconductor.org/packages/3.2/bioc/html/affy.html
# Paquete que contiene funciones para el análisis de datos de microarrays de oligonucleótidos (Affymetrix)
if (!("affy" %in% installed.packages())) { 
  source("http://www.bioconductor.org/biocLite.R"); 
  biocLite("affy");
}


## AFFYPLM
# http://www.bioconductor.org/packages/release/bioc/html/affyPLM.html
# Paquete que extiende y mejora la funcionalidad del paquete affy (herramientas de análisis de calidad)
if (!("affyPLM" %in% installed.packages())) { 
  source("http://www.bioconductor.org/biocLite.R"); 
  biocLite("affyPLM");
}


## GENEFILTER
# http://www.bioconductor.org/packages/release/bioc/html/genefilter.html
# Paquete que proporciona funciones para el filtrado de genes a partir de un conjunto de microarrays
if (!("genefilter" %in% installed.packages())) { 
  source("http://www.bioconductor.org/biocLite.R"); 
  biocLite("genefilter");
}


# LIMMA
# http://www.bioconductor.org/packages/3.2/bioc/html/limma.html
# Paquete para el análisis de expresión diferencial en microarrays, aunque contiene mucha funcionalidad adicional para el análisis de microarrays
if (!("limma" %in% installed.packages())) { 
  source("http://www.bioconductor.org/biocLite.R"); 
  biocLite("limma");
}


# hgu133plus2
# http://www.bioconductor.org/packages/3.2/data/annotation/html/hgu133plus2.db.html
# Paquete de anotación -> nos indica la correspondencia entre las sondas y el gen al que representa, entre otro tipo de información de anotación: términos GO, Entrez IDs, cromosoma, etc
if (!("hgu133plus2.db" %in% installed.packages())) { 
  source("http://www.bioconductor.org/biocLite.R"); 
  biocLite("hgu133plus2.db");
}

# KEGGREST
# http://www.bioconductor.org/packages/release/bioc/html/KEGGREST.html
# Paquete que proporciona una interfaz con el servidor KEGG REST
if (!("KEGGREST" %in% installed.packages())) { 
  source("http://www.bioconductor.org/biocLite.R"); 
  biocLite("KEGGREST");
}

## ------------------------------------------------------------------------

##### 2 . Cargamos los paquetes instalados previamente y ficheros auxiliares
library(GEOquery)
library(affy)
library(affyPLM)
library(genefilter)
library(limma)
library(hgu133plus2.db)
library(KEGGREST)



targets <- readTargets("list.txt", row.names="FileName")


## ------------------------------------------------------------------------

#### 5. Lectura de microarrays (ficheros CEL) y su informacion fenotipica mediante el paquete affy
celfiles <- ReadAffy(filenames=targets$FileName)
# celfiles es un objeto de tipo AffyBatch
class(celfiles)
# Al visualizar el contenido del objeto o variable celfiles, vemos que nos indica que el conjunto tiene 12 muestras o microarrays de Affymetrix
# del tipo hgu133plus2, donde estan descritos 54675 probesets. El tama?o de los arrays es de 1164*1164. La informacion de anotacion es la descrita por hgu133plus2
celfiles

# Acceso a los nombres de las muestras o microarrays
sampleNames(celfiles)
# Acceso a la informacionn fenotipica
pData(celfiles)

# Acceso a los valores de expresion de las sondas
eset <- exprs(celfiles)
dim(eset)
head(eset) # Echamos un vistazo a los valores de intesidad de algunas sondas
# La matriz tiene una dimensionn de 1354896 filas y 12 columnas. Cada fila representa una sonda (probe) y cada columna un microarray determinados


## ------------------------------------------------------------------------
#### 6. Analisis de calidad de los datos crudos
# Veamos algunas de las herramientas: pseudo-imagenes, boxplots e histogramas

## Pseudo-imagenes
# Nos muestra la distribucion espacial de los datos de todos los microarrays del conjunto de datos. 
# Como se puede observar, no se observan artefactos espaciales
image(celfiles)

## Boxplots
# Como ser puede apreciar, los chips tienen una distribucion de valores de intensidad similar
boxplot(celfiles, las=2)

## Histogramas
# Como podemos observar, en el experimento vemos distribuciones similares (aunque con ligeras diferencias) entre los diferentes chips que componen el experimento
hist(celfiles)


## ------------------------------------------------------------------------
#### 7. Pre-procesamiento.
# Objetivo: equilibrar los niveles de intensidad de las muestras de todo el experimento a la vez que se mantiene 
# el efecto debido al tratamiento bajo investigacion
# Etapas (en el caso de microarrays de Affymetrix tipo 3'IVT como el del ejercicio): correccion de ruido de fondo,
# normalizacion, correccion PM-MM y agregacion
# La funcionn rma del paquete affy realiza todas las etapas: 
  # Correccion de background: metodo rma
  # Normalizacion: metodo quantile
  # Correccion PM-MM: PMonly
  # Agregacion: metodo median-polish
eset <- expresso(celfiles,
                 bg.correct = TRUE, 
                 bgcorrect.method="rma",
                 normalize = TRUE, 
                 normalize.method="quantiles", 
                 pmcorrect.method="pmonly", 
                 summary.method="medianpolish",
                 verbose = TRUE) 

# Accedemos a los datos de expresion de los microarrays y comprobamos su dimensin
exprs(eset)
dim(eset)

# La matrix tiene una dimension de 54675 filas y 12 columas. Cada fila representa una probeset y cada columna un microarray determinados
# Vemos como se ha hecho efectiva la etapa de agregacion: los valores de intensidad de todas las sondas que forman parte de un transcrito, deben ser tenidos en cuenta ("agregados") para definir el valor de expresion del gen/probeset.
# Los valores de intensidad de los datos pre-procesados est?n en escala log2 como podemos ver a continuacion:
head(eset)

# Ahora comprobaremos con las  herramientas de analisis de calidad (boxplot e histogramas) indicadas arriba, que los datos pre-procesados tienen mejor
# calidad que los datos crudos
exprseset <- as.data.frame(exprs(eset))		
boxplot(data.frame(exprseset),
        main="Boxplot After Normalization (log scale)",
        col = "white", las = 2)

## ------------------------------------------------------------------------
#### 8. Analisis de expresion diferencial
# Objetivo: identificar aquellos genes cuyos patrones de expresionn difieren, de forma significativa, de acuerdo al fenotipo o a las condiciones experimentales. 
# Lo vamos a realizar mediante el uso de tests estadisticos 'moderados' haciendo uso del paquete limma. Esta etapa esta compuesta de numerosos pasos que detallamos a continuacin

# En primer lugar, vamos a realizar un filtrado no especifico tal y como explicamos en los materiales, previo al analisis de expresion diferencial
# Este filtrado no especifico consiste en eliminar aquellos genes del chip con muy poca variabilidad de expresion entre condiciones. 
# Objetivo:  (1) eliminar genes que con toda probabilidad no van a ser la respuesta a nuestra hipotesis experimental y 
#            (2) reducir el numero de tests a realizar en el analisis de expresion diferencial.

# Sobre este filtro no especifico:
#  - Eliminamos las probesets de control del chip que no son de interes para el investigador: comienzan por AFFX (feature.exclude="^AFFX")
#  - Como medida  de dispersionn para el filtrado por varianza utilizamos IQR (interquartile range) (var.func=IQR) y elegimos el cuartil 0.5 de los valores IQR como umbral (var.cutoff=0.5). Aunque podria aplicarse un filtro menos agresivo
#  - Exigimos ademas que no se eliminen probesets que no tengan identificador Entrez Gene ID (pueden ayudar a construir el modelo estadistico de expresion diferencial si cumplen los valores anteriores)
# nsFilter devuelve una lista. El nuevo objeto ExpressionSet obtenido tras el filtrado est? accesible a trav?s del elemento 'eset'


# La matriz tiene una dimension de 10094 filas y 12 columas. ?si, hemos eliminado 54675-10094 = 44581 probesets



# Recordemos la informacion fenotipica
phenodata
samples <- targets$Classes
samples
# Guardamos la informacionn fenotipica como un tipo de datos 'factor'
samples <- as.factor(samples)

# Creamos la matriz de dise?o teniendo en cuenta las diferentes muestras. Cada fila representa un array y el valor 1 en la columna 
#  nos indica a que tipo de muestra pertenece el array
design <- model.matrix(~0+samples) 
design
colnames(design)
# Proporcionamos a la cabecera de las columnas nombres mas intuitivos 
colnames(design) <- c("hpb_dmso", "hpb_sahm", "kopt_dmso", "kopr_sahm")
design


# Dado un conjunto de arrays pre-procesados y filtrados y el dise?o, ajustamos un modelo lineal a cada gen
# en nuestro conjunto de arrays
fit = lmFit(eset, design)


# Imaginemos que queremos estudiar, unicamente, la expresion diferencial entre arrays del coroides y arrays de la retina (coroides vs retina). 
# Para ello, dise?amos la siguiente matriz de contrastes:
contrast.matrix = makeContrasts(
  hpb_dmsovshpb_sahm = hpb_dmso - hpb_sahm, 
  levels = design)
contrast.matrix

# Podriamos a?dir a la anterior matriz de contrastes tantas comparativas o contrastes como deseemos, por ejemplo, si estamos interesados en las comparativas 
# choroid vs retina, choroid vs huvec y choroid vs iris, nuestra matriz de los 3 contrastes seria la indicada a continuacion: 
#contrast.matrix = makeContrasts(
#  choroid_retina = choroid - retina, 
#  choroid_huvec = choroid - huvec, 
#  choroid_iris <- choroid - iris, 
#  levels = design)
#


# Calculamos una estimacion de los coeficientes y los errores para un conjunto de contrastes
fit.cont <- contrasts.fit(fit, contrast.matrix)
# Calculo del estadistico t moderado y probabilidades de expresion diferencial siguiendo un modelo Bayesiano
res.limma <- eBayes(fit.cont)

# Generacion de la lista de los genes/probesets mas expresados diferencialmente ordenadas por p-value ('sort.by=p'). Se realiza un 
# ajuste de los p-values por el metodo FDR (Benjamini Hochberg) (adjust.method='BH'). El maximo numero de probesets a mostrar 
# es el contenido en nuestro conjunto de arrays (number=nrow(res.limma)) y unicamente se muestran aquellos genes cuyo valor 
# valor p ajustado es inferior a 0.05 y, a la vez, tenga un log fold-change en valor absoluto de al menos 1 (lfc=1) que 
# corresponde genes que est?n sobre-expresados o inhibidos por al menos un factor de 2
# La tabla de resultados muestra la siguiente informacion:
# La primera columna corresponde al identificador del probeset, suele tener formato del tipo <valor_numerico>_at. A continuacion, se muestra 
# la siguiente informacion para cada probeset/gen (aunque un gen puede estar en mas de un probeset)
#  logFC: valor de fold change en base logaritmica (log2)
#  AveExpr: valor de expresion medio entre todos los microarrays (en escala log2)
#  t: estadistico t moderado
#  p.value: valor p asociado al test estadistico
#  adj.P.value: valor p ajustado como resultado de la aplicacion del ajuste de test multiple
#  B: probabilidad en base logaritmica de que el gen est? diferencialmente expresado (cuanto mayor, mejor)
# Ver manual de Limma para conocer mas detalles


results <- topTable(res.limma, p.value=0.05, adjust.method="BH", sort.by="p",number=nrow(res.limma),lfc=1)
# El argumento de topTable coef=<valor> nos indica qu? contraste estudiar: si <valor>=1, devuelve resultados del primer contraste de la matriz de contrastes, si <valor>=2, devuelve resultados del segundo contraste de la matriz de contrastes, etc
# Si solo hay un contraste en nuestra matriz de contrastes, no hace falta indicar ningun valor para coef

dim(results)
# Tenemos un total de 164 probesets con un valor p-ajustado <0.05 y con al menos un lfc>= 1 (valor absoluto) (genes sobre-expresados / inhibidos por un factor de al menos 2)

head(results)



## ------------------------------------------------------------------------
#### 9. An?lisis de alto nivel
# En esta ultima etapa, podemos realizar analisis adicionales que traten de dar una respuesta biologica a nuestro
# conjunto de probesets/genes expresados diferencialmente:
#   - Anotacion
#   - Analisis de agrupamiento
#   - Analisis de enriquecimiento

# Nos vamos a detener en la parte de anotacion, donde vamos a dotar de significado biologico a los genes/probesets que han sido
# identificados como diferencialmente expresados. Por ejemplo, podemos asignar al probeset informacion sobre:
#  Genesymbol
#  EntrezID
#  Pathway
#  Termino GO (Gene Ontology)
# .....

# Veamos las primeras filas de results
head(results)
# Los identificadores de las probesets corresponde a los nombres de filas del data.frame
probesets<-rownames(results)

# Vamos a obtener, para todos los identificadores de probesets del data.frame de resultados, el simbolo del 
# gen al que representa, haciendo uso de la anotacionn del chip hgu133plus2
symbols_list <- mget(probesets,hgu133plus2SYMBOL)
# Pasamos de tipo lista a tipo caracter
symbols<-unlist(symbols_list)
# Vemos a continuacion el SYMBOL para algunos probesets
head(symbols)
# Lo a?adimos a la tabla de resultados
results <- cbind(results, symbols)
# Vemos las primeras filas de results. Se ha a?adido una columna de symbols

head(results)

# Vamos a obtener, para todos los identificadores de probesets del data.frame de resultados, la descripcion del 
# gen al que representa, haciendo uso de la anotacionn del chip hgu133plus2
genename_list <-mget(probesets,hgu133plus2GENENAME)
genename<-unlist(genename_list)
# Vemos a continuacion el genename para algunos probesets
head(genename)
# Lo a?adimos a la tabla de resultados
results <- cbind(results, genename)
head(results)

# Vamos a obtener, para todos los identificadores de probesets del data.frame de resultados, el/los pathway(s) del 
# gen al que representa, haciendo uso de la anotacion del chip hgu133plus2

pathway_list <- mget(probesets,hgu133plus2PATH)
head(pathway_list)


#Funcion para obtener informacion de pathways
get_info_kegg<-function(pathway_ids)
{
  
  info_kegg<-c()
  
  if (any(is.na(pathway_ids))==FALSE) # Si existen identificadores de pathways
  {
    for(current_pathway in pathway_ids) # Recorremos todos los pathway ids
    {
      print(paste('Getting pathway information for',current_pathway))
      # Funci√≥n de KEGGREST que devuelve para un pathway ID, su descripci√≥n (atributo NAME):
      kegg_name<- keggGet(paste('hsa',current_pathway,sep=''))[[1]]$NAME 
      # Concatenamos el resultado en la forma <pathway_id>:info_pathway
      info_kegg<-cbind(info_kegg,paste(current_pathway,kegg_name,sep=':'))
    }
  }
  
  return(paste(info_kegg,collapse=' // '))
  
}
# Obtenemos informacion de los pathways por medio de la funcion get_info_kegg 
kegg_full<-sapply(pathway_list,get_info_kegg)
# Lo a?adimos a la tabla de resultados
results <- cbind(results, kegg_full)
head(results)


# Guardamos los resultados en un fichero que, posteriormente, podemos abrir con un programa de hoja de calculo, por ejemplo Excel o Calc
write.table(results,file='comparison.txt',sep='\t',quote=FALSE,col.names=NA)

# Tambien podemos realizar una tarea similar con el objetvo hgu133plus2GO para obtener las anotaciones de Gene Ontology de un conjunto de probesets

# Ver todas las posibilidades de anotacion de hgu133plus en 
#   http://www.bioconductor.org/packages/release/data/annotation/manuals/hgu133plus2.db/man/hgu133plus2.db.pdf


