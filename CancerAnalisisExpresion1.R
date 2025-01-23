# Importamos las librerias requeridas para el análisis de expresión
library(DESeq2)
library(pheatmap)
library(EnhancedVolcano)
library(ggplot2)

# Importamos la base de datos que se va a analizar

DatosCrudos<- read.csv("DataAnalisis.csv",sep="\t")
head(DatosCrudos)

# Filtramos las columnas de la base de datos de nuestro interés


Data <- DatosCrudos[, c("mRNA", "Control_1", "Control_2", "Control_3", "DCC.2036_1", "DCC.2036_2", "DCC.2036_3")]

countData <- Data[,-1]

# Establecemos las condiciones experimentales para el análisis de expresión
Condiciones<- factor(c("Control","Control","Control","Experimental",
                       "Experimental","Experimental"))

# Convertimos los datos crudos en números enteros (obligatorio)
countData <- round(countData)

# Convertimos en un Objeto DESeq2 para que pueda ser analizado
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = data.frame(Condiciones),
                              design = ~ Condiciones)

# Con esta función realizamos el análisis de expresión (incluye un pipeline)
dds <- DESeq(dds)

# Extraemos los resultados en una tabla del pipeline
res <- results(dds)
head(res)

# Podemos realizar algunos análisis de los resultados obtenidos por filtración 

# Antes convertimos los resultados extraídos en un objeto dataframe

DataframeDatos<- as.data.frame(res)

# Le agregamos la columna con el nombre de cada transcrito (el orden sigue igual)
DataframeDatos$Transcritos<- DatosCrudos$mRNA

DataframeDatos <- DataframeDatos[, c("Transcritos", setdiff(colnames(DataframeDatos), "Transcritos"))]

# Realizamos un procesamiento, eliminando todos los valores NA del dataframe 
DataframeLimpio <- na.omit(DataframeDatos)

# Filtramos qué transcritos presentaron un valor de log2FoldChange mayor a 6

Transcritos_Mayores_2<- DataframeLimpio[abs(DataframeLimpio$log2FoldChange)>6,]

Transcritos_Mayores_2

# Representar visualmente los resultados de Análisis de expresión diferencial

# Realizar algunos ajustes para reducir el ruido biológico
resLFC <- lfcShrink(dds, coef="Condiciones_Experimental_vs_Control", type="apeglm")
resNormal<- lfcShrink(dds, coef="Condiciones_Experimental_vs_Control",type="normal")

# PlotMA
plotMA(res, ylim=c(-4,4))
plotMA(resLFC, ylim=c(-4,4))
plotMA(resNormal, ylim=c(-4,4))

# VolcanoPlot para genes altamente significativos con alta variación en FoldChange
EnhancedVolcano(res,
                lab = DatosCrudos$mRNA,  # Etiquetas de genes
                x = 'log2FoldChange',
                y = 'padj',       
                pCutoff = 0.05,       # Umbral de p-valor
                FCcutoff = 1,         # Umbral de cambio de pliegue
                title = 'Volcano Plot',
                caption = 'Análisis de DESeq2',
                col = c('gray', 'red', 'blue', 'green'),  # Colores para los puntos
                legendPosition = 'top',
                legendLabSize = 12,
                legendIconSize = 4)

# Creación de un gráfico de PCA (Análisis de componentes principales)

rld <- rlog(dds, blind = TRUE)

pca_plot <- plotPCA(rld, intgroup = "Condiciones", returnData = TRUE)
percentVar <- round(100 * attr(pca_plot, "percentVar"))

ggplot(pca_plot, aes(x = PC1, y = PC2, color = Condiciones)) +
  geom_point(size = 4) +  
  labs(
    title = "PCA: Análisis de Componentes Principales",
    x = paste0("PC1 (", percentVar[1], "% varianza)"),
    y = paste0("PC2 (", percentVar[2], "% varianza)"),
    color = "Condiciones"
  ) +
  theme_minimal(base_size = 15) + 
  theme(
    plot.title = element_text(hjust = 0.5), 
    legend.position = "top" 
  )  

# Creación de un mapa de calor 

# A los resultados de expresión le agregamos la columna de nombres de transcritos
res$mRNA <- DatosCrudos$mRNA[1:nrow(res)]

# Filtramos y eliminamos los valores NAs
res<- na.omit(res)

# Seleccionamos los primeros 50 genes más significativos
res_sig <- res[res$padj < 0.05, ][1:50,]
head(res_sig)

# Nos vamos a los datos crudos y seleccionamos esos 50 genes
genes_seleccionados <- DatosCrudos[DatosCrudos$mRNA%in% res_sig$mRNA, ]

# Establecemos las variables a analizar en el mapa de calor
datos_para_heatmap <- genes_seleccionados[, c("Control_1", "Control_2", "Control_3", 
                                              "DCC.2036_1", "DCC.2036_2", "DCC.2036_3")]

# Generamos el mapa de calor 
pheatmap(datos_para_heatmap, 
         cluster_rows = TRUE, # Agrupar las filas (genes)
         cluster_cols = TRUE, # Agrupar las columnas (condiciones)
         show_rownames = FALSE, # No mostrar los nombres de las filas para mayor claridad
         show_colnames = TRUE, # Mostrar los nombres de las condiciones
         scale = "row", # Escalar por filas (genes) para que las diferencias entre los genes sean más evidentes
         color = colorRampPalette(c("blue", "white", "red"))(50), # Colores del mapa de calor
         fontsize = 10) # Ajustar el tamaño de la fuente
