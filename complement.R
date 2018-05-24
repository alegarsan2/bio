
###### 1. InstalaciÃ³n de paquetes necesarios para la realizaciÃ³n de este ejercicio

# Geo Query
# http://www.bioconductor.org/packages/release/bioc/html/GEOquery.html
# Paquete que permite descargar del repositorio GEO (Gene Expression Omnibus) conjuntos de microarrays
if (!("GEOquery" %in% installed.packages())) { 
  source("http://www.bioconductor.org/biocLite.R"); 
  biocLite("GEOquery");
}


# AFFY
# http://www.bioconductor.org/packages/3.2/bioc/html/affy.html
# Paquete que contiene funciones para el anÃ¡lisis de datos de microarrays de oligonucleÃ³tidos (Affymetrix)
if (!("affy" %in% installed.packages())) { 
  source("http://www.bioconductor.org/biocLite.R"); 
  biocLite("affy");
}


## AFFYPLM
# http://www.bioconductor.org/packages/release/bioc/html/affyPLM.html
# Paquete que extiende y mejora la funcionalidad del paquete affy (herramientas de anÃ¡lisis de calidad)
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
# Paquete para el anÃ¡lisis de expresiÃ³n diferencial en microarrays, aunque contiene mucha funcionalidad adicional para el anÃ¡lisis de microarrays
if (!("limma" %in% installed.packages())) { 
  source("http://www.bioconductor.org/biocLite.R"); 
  biocLite("limma");
}


# hgu133plus2
# http://www.bioconductor.org/packages/3.2/data/annotation/html/hgu133plus2.db.html
# Paquete de anotaciÃ³n -> nos indica la correspondencia entre las sondas y el gen al que representa, entre otro tipo de informaciÃ³n de anotaciÃ³n: tÃ©rminos GO, Entrez IDs, cromosoma, etc
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

# Cargar fichero auxiliar 'get_info_kegg.R'. Sustituye la línea de abajo por el directorio donde almacenes este fichero
source('/path_to/get_info_kegg.R') # Funcion que obtiene informacio de pathways a partir de una lista de identificadores
# Si utilizas Windows, debes poner source('C:/path_to/get_info_kegg.R')


## ------------------------------------------------------------------------

#### 3. Descarga del conjunto de datos haciendo uso de GEOquery

## El conjunto de datos a analizar tiene como referencia en GEO GSE20986 (http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE20986). 
## El trabajo publicado con este conjunto de datos se titula "Comparative gene expression profiling of human umbilical vein endothelial cells and ocular vascular endothelial cells"
# El objetivo del trabajo es investigar las diferencias entre células endoteliales de vena umbilical humana y las células endoteliales vasculares oculares (coroides, retina e iris) con el 
# propósito de determinar si estas diferencias pueden mejorar el entendimiento de enfermedades oculares

# Descargamos el conjunto de datos, creando un directorio llamado GSE20986 contenido en el directorio actual o de trabajo. Dentro del directorio GSE20986
# encontramos el fichero GSE20986_RAW.tar que contiene los 12 ficheros CEL (tambiÃ©n comprimidos) y un fichero filelist.txt con el listado de los
# ficheros descargados
x = getGEOSuppFiles("GSE20986")
x

## Descomprimir el fichero descargado en un directorio llamado "data" contenido en el directorio de trabajo
untar("GSE20986/GSE20986_RAW.tar", exdir = "data")
## Si nos fijamos dentro de "data", todos los ficheros CEL estan tambien comprimidos

## Listar todos los ficheros contenidos en el directorio data con el sufijo 'gz' (fichero comprimido) y guardarlos en la variable cels
cels = list.files("data/", pattern = "[gz]")
## Mediante sapply, aplicamos la descompresion (gunzip) a todos los ficheros comprimidos contenidos en ./data, obteniendo todos los ficheros .CEL
sapply(paste("data", cels, sep = "/"), gunzip)

## ------------------------------------------------------------------------

#### 4. Construir la informacionn fenotipica del conjunto de datos. Tenemos 12 microarrays:
  # 3 microarrays correspondientes a muestras de celulas endoteliales de la vena umbilical (replicas biologicas) -> huvec
  # 3 microarrays correspondientes a muestras de celulas endoteliales vasculares del coroides (replicas biologicas) -> choroid
  # 3 microarrays correspondientes a muestras de celulas endoteliales vasculares de la retina (replicas biologicas) -> retina
  # 3 microarrays correspondientes a muestras de celulas endoteliales vasculares del iris (replicas biologicas) -> iris
phenodata = matrix(rep(list.files("data"), 2), ncol =2)
class(phenodata)
phenodata <- as.data.frame(phenodata)
colnames(phenodata) <- c("Name", "FileName")
phenodata$Targets <- c("iris", 
                       "retina", 
                       "retina", 
                       "iris", 
                       "retina", 
                       "iris", 
                       "choroid", 
                       "choroid", 
                       "choroid", 
                       "huvec", 
                       "huvec", 
                       "huvec")

# La informacionn fenotipica nos muestra la correspondencia entre el fichero CEL y el origen de las celulas (cordonn umbilical, coroides, retina o iris)
# Guardamos la informacionn fenotpica en el fichero phenodata.txt contenido en el directorio data
write.table(phenodata, "phenodata.txt", quote = F, sep = "\t", row.names = F)


## ------------------------------------------------------------------------

#### 5. Lectura de microarrays (ficheros CEL) y su informacion fenotipica mediante el paquete affy
celfiles <- ReadAffy(phenoData = "phenodata.txt", celfile.path="data")
# celfiles es un objeto de tipo AffyBatch
class(celfiles)
# Al visualizar el contenido del objeto o variable celfiles, vemos que nos indica que el conjunto tiene 12 muestras o microarrays de Affymetrix
# del tipo hgu133plus2, donde estan descritos 54675 probesets. El tamaño de los arrays es de 1164*1164. La informacion de anotacion es la descrita por hgu133plus2
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
# A partir de este punto, tened a mano los materiales de la asignatura, seccion "2.2. Flujo de analisis de datos de expresion genica mediante microarrays"
# Como recordareis, un flujo de analisis tipico de expresión diferencial de microarrays consiste en: analisis de calidad, pre-procesamiento, analisis de expresion diferencial y analisis de alto nivel

## ------------------------------------------------------------------------
#### 6. Analisis de calidad de los datos crudos
# Veamos algunas de las herramientas: pseudo-imagenes, boxplots e histogramas

## Pseudo-imagenes
# Nos muestra la distribucion espacial de los datos de todos los microarrays del conjunto de datos. 
# Cada vez que pulsemos return nos mostrara un chip del conjunto hasta mostrarlos todos
# Como se puede observar, no se observan artefactos espaciales
image(celfiles)

## Boxplots
# Los boxplots o diagramas de cajas, son unas herramientas que proporcionan de forma visual la distribucion de los datos 
# (valores de intensidad de las sondas en este caso) de cada microarray del experimento a traves de cinco medidas: valor minimo,
# cuartil Q1, mediana o cuartil Q2, cuartil Q3 y valor maximo
# Como ser puede apreciar, los chips tienenr una distribucion de valores de intensidad similar
boxplot(celfiles, las=2)

## Histogramas
# los histogramas se utilizan para obtener de forma grafica la distribucion de los valores de intensidad de las sondas de cada uno de los 
# microarrays del experimento. Es una herramienta complementaria a los boxplots que nos proporciona otra vision de la distribucion de los valores de intensidad de las sondas
# Como podemos observar, en el experimento vemos distribuciones similares (aunque con ligeras diferencias) entre los diferentes chips que componen el experimento
hist(celfiles)

## MAplots
# Nos los vamos a ver en este ejercicio guiado, pero los podeis obtener mediante la funcin MAplot del paquete AffyPLM (http://www.bioconductor.org/packages/3.2/bioc/html/affyPLM.html)

# POr tanto, de acuerdo a estas herramientas de calidad, todos los arrays, en general, tienen buena calidad.

## Para herramietnas mas complejas de analisis de calidad, echad un vistazo a la siguiente bibliografia:
# Bioconductor Case Studies: http://www-huber.embl.de/pub/pdf/HahneHuberGentlemanFalcon2008.pdf


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
celfiles.rma <- rma(celfiles)
celfiles.rma

# Los datos de microarrays pre-procesados se encuentran en la variable celfiles.rma, que es un objeto tipo ExpressionSet
class(celfiles.rma)

# Accedemos a los datos de expresion de los microarrays y comprobamos su dimensin
eset<-exprs(celfiles.rma)
dim(eset)
# La matrix tiene una dimension de 54675 filas y 12 columas. Cada fila representa una probeset y cada columna un microarray determinados
# Vemos como se ha hecho efectiva la etapa de agregacion: los valores de intensidad de todas las sondas que forman parte de un transcrito, deben ser tenidos en cuenta ("agregados") para definir el valor de expresion del gen/probeset.
# Los valores de intensidad de los datos pre-procesados están en escala log2 como podemos ver a continuacion:
head(eset)

# Ahora comprobaremos con las  herramientas de analisis de calidad (boxplot e histogramas) indicadas arriba, que los datos pre-procesados tienen mejor
# calidad que los datos crudos
boxplot(celfiles.rma,las=2)
hist(celfiles.rma)

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
# nsFilter devuelve una lista. El nuevo objeto ExpressionSet obtenido tras el filtrado está accesible a través del elemento 'eset'

celfiles.rma_filtered<-nsFilter(celfiles.rma, require.entrez=FALSE, var.func=IQR, var.cutoff=0.5, feature.exclude="^AFFX")$eset
exprdata_filtered<-exprs(celfiles.rma_filtered)
dim(exprdata_filtered)
# La matriz tiene una dimension de 10273 filas y 12 columas. Âsi, hemos eliminado 54675-10273 = 44402 probesets



# Recordemos la informacion fenotipica
phenodata
samples <- celfiles$Targets
samples
# Guardamos la informacionn fenotipica como un tipo de datos 'factor'
samples <- as.factor(samples)

# Creamos la matriz de diseño teniendo en cuenta las diferentes muestras. Cada fila representa un array y el valor 1 en la columna 
#  nos indica a que tipo de muestra pertenece el array
design <- model.matrix(~0+samples) 
design
colnames(design)
# Proporcionamos a la cabecera de las columnas nombres mas intuitivos 
colnames(design) <- c("choroid", "huvec", "iris", "retina")
design


# Dado un conjunto de arrays pre-procesados y filtrados y el diseño, ajustamos un modelo lineal a cada gen
# en nuestro conjunto de arrays
fit = lmFit(celfiles.rma_filtered, design)


# Imaginemos que queremos estudiar, unicamente, la expresion diferencial entre arrays del coroides y arrays de la retina (coroides vs retina). 
# Para ello, diseñamos la siguiente matriz de contrastes:
contrast.matrix = makeContrasts(
  choroid_retina = choroid - retina, 
  levels = design)
contrast.matrix

# Podriamos añdir a la anterior matriz de contrastes tantas comparativas o contrastes como deseemos, por ejemplo, si estamos interesados en las comparativas 
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
# corresponde genes que estén sobre-expresados o inhibidos por al menos un factor de 2
# La tabla de resultados muestra la siguiente informacion:
# La primera columna corresponde al identificador del probeset, suele tener formato del tipo <valor_numerico>_at. A continuacion, se muestra 
# la siguiente informacion para cada probeset/gen (aunque un gen puede estar en mas de un probeset)
#  logFC: valor de fold change en base logaritmica (log2)
#  AveExpr: valor de expresion medio entre todos los microarrays (en escala log2)
#  t: estadistico t moderado
#  p.value: valor p asociado al test estadistico
#  adj.P.value: valor p ajustado como resultado de la aplicacion del ajuste de test multiple
#  B: probabilidad en base logaritmica de que el gen esté diferencialmente expresado (cuanto mayor, mejor)
# Ver manual de Limma para conocer mas detalles


results <- topTable(res.limma, p.value=0.05, adjust.method="BH", sort.by="p",number=nrow(res.limma),lfc=1)
# El argumento de topTable coef=<valor> nos indica qué contraste estudiar: si <valor>=1, devuelve resultados del primer contraste de la matriz de contrastes, si <valor>=2, devuelve resultados del segundo contraste de la matriz de contrastes, etc
# Si solo hay un contraste en nuestra matriz de contrastes, no hace falta indicar ningun valor para coef

dim(results)
# Tenemos un total de 1701 probesets con un valor p-ajustado <0.05 y con al menos un lfc>= 1 (valor absoluto) (genes sobre-expresados / inhibidos por un factor de al menos 2)

# Vamos a fijarnos en las primeras filas de esta tabla
head(results)
# En la matriz de contrastes, indicamos choroid - retina, es decir, choroid vs retina. Por tanto aquellos probesets cuyo logFC o estadistico t
# tenga un valor positivo, significa sobre-expresion de ese probeset en choroid con respecto a retina. Por el contrario, un valor 
# negativo significa inhibicionn de ese probeset en choroid con respecto a retina

# Ayuda: En el manual completo de LIMMA, podemos obtener mas informacionn (http://www.bioconductor.org/packages/release/bioc/vignettes/limma/inst/doc/usersguide.pdf)
# El capitulo 8.2 nos muestra un ejemplo sencillo con arrays de Affymetrix



## ------------------------------------------------------------------------
#### 9. Análisis de alto nivel
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
# Lo añadimos a la tabla de resultados
results <- cbind(results, symbols)
# Vemos las primeras filas de results. Se ha añadido una columna de symbols

head(results)

# Vamos a obtener, para todos los identificadores de probesets del data.frame de resultados, la descripcion del 
# gen al que representa, haciendo uso de la anotacionn del chip hgu133plus2
genename_list <-mget(probesets,hgu133plus2GENENAME)
genename<-unlist(genename_list)
# Vemos a continuacion el genename para algunos probesets
head(genename)
# Lo añadimos a la tabla de resultados
results <- cbind(results, genename)
head(results)

# Vamos a obtener, para todos los identificadores de probesets del data.frame de resultados, el/los pathway(s) del 
# gen al que representa, haciendo uso de la anotacion del chip hgu133plus2

pathway_list <- mget(probesets,hgu133plus2PATH)
head(pathway_list)

# Obtenemos informacion de los pathways por medio de la funcion get_info_kegg (get_info_kegg.R)
kegg_full<-sapply(pathway_list,get_info_kegg)
# Lo añadimos a la tabla de resultados
results <- cbind(results, kegg_full)
head(results)


# Guardamos los resultados en un fichero que, posteriormente, podemos abrir con un programa de hoja de calculo, por ejemplo Excel o Calc
write.table(results,file='choroidal_vs_retinal_endothelial_cells.txt',sep='\t',quote=FALSE,col.names=NA)

# Tambien podemos realizar una tarea similar con el objetvo hgu133plus2GO para obtener las anotaciones de Gene Ontology de un conjunto de probesets

# Ver todas las posibilidades de anotacion de hgu133plus en 
#   http://www.bioconductor.org/packages/release/data/annotation/manuals/hgu133plus2.db/man/hgu133plus2.db.pdf


## ------------------------------------------------------------------------
#### Consideraciones finales
### Este ejercicio guiado ha sido diseñado con propositos docentes. Los microarrays utilizados en este ejemplo (hgu133plus2), se siguen
# utilizando en la actualidad (son del tipo 3'IVT). Sin embargo, existen otros microarrays de Affymetrix, tal y como vimos en los materiales, mucho mas sofisticados
# y mas completos, como los whole-transcript. Asi, encontramos el modelo Human GEne ST 2.0 array. Existen paquetes en bioconductor para el procesamiento de estos arrays, como
# por ejemplo el paquete oligo (http://www.bioconductor.org/packages/release/bioc/html/oligo.html). El workflow de analisis es muy similar al visto
# en este ejercicio guiado.


## ------------------------------------------------------------------------
## ------------------------------------------------------------------------
### Bibliografia
#   - Bioconductor Case Studies: http://www-huber.embl.de/pub/pdf/HahneHuberGentlemanFalcon2008.pdf
#   - Bioinformatics and Computational Biology Solutions using R and Bioconductor: http://link.springer.com/book/10.1007%2F0-387-29362-0
