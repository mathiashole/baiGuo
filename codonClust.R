#Load the required libraries

pkgs <- sort(c('tidyverse', 'factoextra', 'cluster', 'mclust', 
               'kernlab'
))

pkgs_install <- pkgs[!(pkgs %in% installed.packages()[,"Package"])]
if(length(pkgs_install)) install.packages(pkgs_install)

library(factoextra)
library(clv)
library(ggplot2)
#library(pheatmap)#Import data
library(tidyverse) # add_column() function
library(mclust) # Gaussian Mixture Models
library(cluster) # Partitioning Around Medoids (PAM)
library(kernlab) # scpectral clustering


codon_usage <- read.csv("/home/usuario/Data_Rstudio/codon_freq/codon_usage.csv")
adenodf <- codon_usage

#Create advir dataframe containing only virus columns
advir <- adenodf[grep(c("vrl"), adenodf$Kingdom), ]#use a matchpattern object containing "Human" and "adenovirus" strings to extract 
#Human adenovirus data
matchpattern <- c("Human", "adenovirus")#Use another grepl to select rows containing "adenovirus" 
advir<- advir[grepl("Human adenovirus", advir$SpeciesName),]#Make Rownames => Species
rownames(advir) <- advir$SpeciesName


# Convertir las columnas que deben contener números a numéricas
columnas_numericas <- c("UUU", "UUC", "UUA", "UUG", "CUU", "CUC", "CUA", "CUG", 
                                     "AUU", "AUC", "AUA", "AUG", "GUU", "GUC", "GUA", "GUG", 
                                     "GCU", "GCC", "GCA", "GCG", "CCU", "CCC", "CCA", "CCG", 
                                     "UGG", "GGU", "GGC", "GGA", "GGG", "UCU", "UCC", "UCA", 
                                     "UCG", "AGU", "AGC", "ACU", "ACC", "ACA", "ACG", "UAU", 
                                     "UAC", "CAA", "CAG", "AAU", "AAC", "UGU", "UGC", "CAU", 
                                     "CAC", "AAA", "AAG", "CGU", "CGC", "CGA", "CGG", "AGA", 
                                     "AGG", "GAU", "GAC", "GAA", "GAG", "UAA", "UAG", "UGA")

count <- advir[,6:69]
  
count[, columnas_numericas] <- apply(count[, columnas_numericas], 2, as.numeric)

# Verificar la estructura de datos después de la conversión
str(count)

# Aplicar prcomp
pca <- prcomp(count, center = TRUE, scale = TRUE)

#pca = prcomp(dmatrix, center = TRUE, scale = TRUE)
summary(pca)

### PLOT PCA ___________________________________________________________________

fviz_pca_biplot(pca,
                label="var",
                habillage = advir$SpeciesName )

transform = as.data.frame(-pca$x[,1:2])



#_______________________________________________________________________________

fviz_nbclust(transform, kmeans, method = 'wss')
fviz_nbclust(transform, kmeans, method = 'silhouette')
gp_st <- fviz_nbclust(transform, kmeans, method = 'gap_stat')



k = 3

## kmean
# 
# kmeans_transform = kmeans(transform, centers = k, nstart = 50)
# fviz_cluster(kmeans_transform, data = transform, ggtheme = theme_minimal())
# 
# trnasform_clustered <- transform |> add_column(cluster = factor(kmeans_transform$cluster))
# 
# ggplot(trnasform_clustered, aes(x = PC1, y = PC2, color = cluster)) + geom_point() + 
#   stat_ellipse(geom="polygon", aes(fill = cluster),
#                alpha = 0.2,
#                show.legend = FALSE,
#                level = 0.95)
# 
# centroids <- as_tibble(kmeans_transform$centers, rownames = "cluster")
# centroids$cluster <- as.factor(centroids$cluster)
# ggplot(transform, aes(x = PC1, y = PC2, color = cluster)) + 
#   geom_point() +
#   geom_point(data = centroids, aes(x = PC1, y = PC2, color = cluster), shape = 3, size = 10)


## HC 

# dpcamatrix <- dist(transform)
# hc_transform <- hclust(dpcamatrix, method = "complete")
# 
# #fviz_dend(hc_transform, k = k)
# clusters <- cutree(hc_transform, k = k)
# 
# ## OPCIONAL PLOT________________________________________________________________
# 
# cluster_complete <- transform |>
#   add_column(cluster = factor(clusters))
# 
# ggplot(cluster_complete, aes(PC1, PC2, color = cluster)) +
#   geom_point()
# 
# ##______________________________________________________________________________
# 
# fviz_cluster(list(data = transform, cluster = cutree(hc_transform, k = k)),  ggtheme = theme_minimal())


## Spectral clustering
#install.packages("kernlab")
# library(kernlab)
# 
# cluster_spec <- specc(as.matrix(transform), centers = k)
# fviz_cluster(list(data = transform, cluster = cluster_spec), geom = "point")
# ggplot(transform |> 
#          add_column(cluster = factor(cluster_spec)),
#        aes(PC1, PC2, color = cluster)) + 
#   geom_point()


## Partitioning Around Medoids (PAM)
# library(cluster)
# 
# d <- dist(transform)
# str(d)
# 
# pam_result <- pam(d, k = k)
# 
# fviz_cluster(c(pam_result, list(data = transform)), ggtheme = theme_minimal())

## Gaussian Mixture Models
# library(mclust)

# mmclust <- Mclust(transform, G = k)
# 
# fviz_cluster(list(data = transform, cluster = mmclust$classification), ggtheme = theme_minimal())



perform_clustering <- function(transform, k, method) {
  if (method == 'kmeans') {
    kmeans_result <- kmeans(transform, centers = k, nstart = 50)
    fviz_cluster(kmeans_result, data = transform, ggtheme = theme_minimal())
  } else if (method == 'hclust') {
    dpcamatrix <- dist(transform)
    hc_result <- hclust(dpcamatrix, method = "complete")
    fviz_cluster(list(data = transform, cluster = cutree(hc_result, k = k)),  ggtheme = theme_minimal())
  } else if (method == 'specc') {
    cluster_spec <- specc(as.matrix(transform), centers = k)
    fviz_cluster(list(data = transform, cluster = cluster_spec), ggtheme = theme_minimal())
  } else if (method == 'pam') {
    dpcamatrix <- dist(transform)
    pam_result <- pam(d, k = k)
    fviz_cluster(c(pam_result, list(data = transform)), ggtheme = theme_minimal())
  } else if (method == 'mclust') {
    mmclust <- Mclust(transform, G = k)
    fviz_cluster(list(data = transform, cluster = mmclust$classification), ggtheme = theme_minimal())
  } else {
    stop("Método de clustering no reconocido.")
  }
}

perform_clustering(transform, k = k, method = 'kmeans')
perform_clustering(transform, k = k, method = 'hclust')
perform_clustering(transform, k = k, method = 'specc')
perform_clustering(transform, k = k, method = 'pam')
perform_clustering(transform, k = k, method = 'mclust')



################################################################################

# fix classification supervizado

# # Instalar y cargar las librerías necesarias
# #install.packages("caret")
# library(caret)
# #install.packages("class")
# library(class)
# #install.packages("e1071")
# library(e1071)
# #install.packages("randomForest")
# library(randomForest)
# #install.packages("nnet")
# library(nnet)
# #install.packages("rpart")
# library(rpart)
# 
# 
# # Cargar el conjunto de datos iris
# data(iris)
# 
# # Definir variables predictoras (X) y variable de respuesta (y)
# x <- iris[, 1:4]
# y <- as.factor(iris$Species)
# 
# # Crear un dataframe para almacenar las precisiones de diferentes modelos
# acc_r <- matrix(0, nrow = 10, ncol = 6)  # 6 modelos en total
# 
# # Crear una función para calcular la precisión
# calculate_accuracy <- function(predicted, actual) {
#   confusion <- confusionMatrix(predicted, actual)
#   return(confusion$overall["Accuracy"])
# }
# 
# for (i in 1:10) {
#   # Creamos un split aleatorio
#   sample_indices <- createDataPartition(y, p = 0.7, list = FALSE)
#   X_train <- x[sample_indices, ]
#   y_train <- y[sample_indices]
#   X_test <- x[-sample_indices, ]
#   y_test <- y[-sample_indices]
#   
#   # Instanciamos los modelos
#   nn1 <- knn(train = X_train, test = X_test, cl = y_train, k = 1)
#   nn3 <- knn(train = X_train, test = X_test, cl = y_train, k = 3)
#   svc <- svm(X_train, y_train)
#   dt <- rpart(y_train ~ ., data = as.data.frame(cbind(X_train, y_train)))
#   rf <- randomForest(X_train, y_train, ntree = 100)
#   lr <- glm(y_train ~ ., data = as.data.frame(cbind(X_train, y_train)), family = "binomial")
#   ann <- neuralnet(y_train ~ ., data = as.data.frame(cbind(X_train, y_train)), hidden = c(5))
#   
#   # Predecimos sobre test con los 6 modelos
#   yhat_nn1 <- knn(train = X_train, test = X_test, cl = y_train, k = 1)
#   yhat_nn3 <- knn(train = X_train, test = X_test, cl = y_train, k = 3)
#   yhat_svc <- predict(svc, newdata = X_test)
#   yhat_dt <- predict(dt, newdata = X_test)
#   yhat_rf <- predict(rf, newdata = X_test)
#   yhat_lr <- predict(lr, newdata = as.data.frame(X_test), type = "response")
#   yhat_ann <- as.numeric(predict(ann, newdata = as.data.frame(X_test)))
#   
#   # Medimos la accuracy de cada modelo en esta iteración
#   acc_r[i, 1] <- calculate_accuracy(yhat_nn1, y_test)
#   acc_r[i, 2] <- calculate_accuracy(yhat_nn3, y_test)
#   acc_r[i, 3] <- calculate_accuracy(yhat_svc, y_test)
#   acc_r[i, 4] <- calculate_accuracy(yhat_dt, y_test)
#   acc_r[i, 5] <- calculate_accuracy(yhat_rf, y_test)
#   acc_r[i, 6] <- calculate_accuracy(round(yhat_lr), y_test)
# }
# 
# # Boxplot para visualizar las precisiones
# boxplot(acc_r, names = c('1-NN', '3-NN', 'SVM', 'Decision Tree', 'Random Forest', 'Logistic Regression'),
#         ylab = 'Accuracy', main = 'Comparación de precisión de modelos')
