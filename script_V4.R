library(dplyr)
library(ggplot2)
library(FactoMineR)
library(factoextra)
library(DiscriMiner)


#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
# Préparation de la BDD----
#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

df <- read.csv2("data/phoneme.csv", sep= ";")
# Sélection des deux dernières colonnes
last_two_cols <- df[, (ncol(df) - 1):ncol(df)]
# On vire les deux dernières colonnes
df <- df[, -c((ncol(df) - 1):ncol(df))]
# Fusion
df <- cbind(last_two_cols, df)
df$g = as.factor(df$g)
df <- df %>% mutate(across(starts_with("x."), as.numeric))
df <- df[, -which(names(df) == "row.names")]
df_train <- df[grepl("train", df$speaker), ]
df_train  = df_train[,-c(2)]
df_test <- df[grepl("test", df$speaker), ]
df_test  = df_test[,-c(2)]
df  = df[,-c(2)]
phoneme_groups <- split(df, df$g)


#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
# Statistiques Descriptives----
#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
# Fonction pour calculer les statistiques descriptives
descriptive_stats <- function(data) {
  # Sélectionner uniquement les colonnes numériques (supposons que les colonnes "x.1" à "x.256" sont numériques)
  numeric_data <- data[, grepl("^x\\.", colnames(data))]
  
  summary_stats <- data.frame(
    mean = colMeans(numeric_data),
    median = apply(numeric_data, 2, median),
    sd = apply(numeric_data, 2, sd),
    min = apply(numeric_data, 2, min),
    max = apply(numeric_data, 2, max)
  )
  return(summary_stats)
}

# Appliquer la fonction aux groupes de phonèmes
phoneme_stats <- lapply(phoneme_groups, descriptive_stats)
# Convertir la liste en data.frame
phoneme_stats_df <- do.call(rbind, phoneme_stats)
# Ajouter une colonne avec les noms des phonèmes
phoneme_stats_df$phoneme <- rep(names(phoneme_groups), each = nrow(phoneme_stats[[1]]))
# Ajouter une colonne avec les numéros de variable (périodogramme)
phoneme_stats_df$variable <- factor(rep(1:nrow(phoneme_stats[[1]]), length(unique(phoneme_stats_df$phoneme))), levels = 1:nrow(phoneme_stats[[1]]))
# Convertir les données en format long
long_phoneme_stats_df <- reshape2::melt(phoneme_stats_df, id.vars = c("phoneme", "variable"), variable.name = "statistic", value.name = "value")
# Ajouter une colonne avec les numéros de variable (périodogramme)
mean_data <- long_phoneme_stats_df[long_phoneme_stats_df$statistic == "mean", ]
# Extraire les moyennes de long_phoneme_stats_df
mean_data <- long_phoneme_stats_df[long_phoneme_stats_df$statistic == "mean",]
# Créer une colonne "variable" avec les numéros de variable (périodogramme) pour chaque phonème
mean_data$variable <- factor(rep(1:nrow(phoneme_stats[[1]]), length(unique(mean_data$phoneme))), levels = 1:nrow(phoneme_stats[[1]]))

library(scales)
palette <- brewer_pal(palette = "Set2")(5)
# Ajouter une colonne avec les numéros de variable (périodogramme)
ggplot(mean_data, aes(x = variable, y = value, color = phoneme, group = phoneme)) +
  geom_line(size=0.9) +
  scale_color_manual(values = palette) +
  theme_minimal() +
  labs(
    x = "Fréquence",
    y = "Valeur moyenne") +
  theme(legend.position = "bottom",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_line(size = 0.2),
        panel.grid.minor.y = element_line(size = 0.1),
        plot.title = element_text(hjust = 0.5,face = "bold"),
        axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"))

# Proportions phonèmes
proportions <- table(df_train$g)/nrow(df_train)

prop2 <- data.frame(categories = names(proportions), proportions = proportions)

ggplot(prop2, aes(x = categories, y = proportions, fill = categories)) + 
  geom_bar(stat = "identity") +
  xlab("Catégories") +
  ylab("Proportions") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

table(df_train$g/nrow(df_train))


#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
# ACP----
#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

res.pca2 = PCA(df_train, quali.sup = 1, axes=c(1,2), graph=FALSE)

fviz_eig(res.pca2,	choice	=	c("variance",	"eigenvalue"),	
         geom	=	c("bar",	"line"),	barfill	=	"royalblue1",
         barcolor	=	"grey35",	linecolor	=	"black",	
         ncp	=	16,	addlabels	=	TRUE,
         main="Histogramme des valeurs propres")

FactoMineR::plot.PCA(res.pca2,choix="ind",habillage=11,label="quali")

# ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
# AFD----
# ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
data.cr<-scale(df_train[,c(-1)],center=TRUE,scale=TRUE)
data.dist<-dist(data.cr)
data.hca<-hclust(data.dist,method="ward.D2")

par(mfrow=c(1,2)) 
# barplot
barplot(rev(data.hca$height)[1:30],main="Evolution du critére d agrégation - VISA")

# Dendrogramme
plot(data.hca,hang=-1,cex=0.8)
res.cutree <- cutree(data.hca,3)
centers <- aggregate(data.cr, by=list(res.cutree), FUN=mean)
centers <- centers[,-1]
data.kmeans <- kmeans(data.cr, centers=centers, algorithm="MacQueen")
data.kmeans.alea <- kmeans(data.cr, centers=centers, nstart = 50, iter.max = 50, algorithm="MacQueen")
data.tot <- cbind(as.data.frame(df_train), Clas = as.factor(data.kmeans$cluster)) 
ncol(data.tot)
table(data.tot[,c(1,258)])
chisq.test(table(data.tot[,c(1,258)])) # Chi square test

# Sélection des variables explicatives
X <- df_train[, 2:257]

# Variable Y
y <- df_train[, 1]

res.desDA <- desDA(X,y, covar = "within")
head(round(res.desDA$power,3),20)
corRatio(res.desDA$scores[,1],y)
FRatio(res.desDA$scores[,1],y)
res.geoDA <- geoDA(X,y,validation="crossval")
res.geoDA$error_rate 
res.geoDA$confusion

pred.app <- classify(res.geoDA,X)
table.BC <- table(y,pred.app$pred_class)
table.BC # matrice de confusion
err.rate.app <- 1- sum(diag(table.BC))/sum(table.BC)
err.rate.app # taux d'erreur
table.BC[2,2]/sum(table.BC[2,]) # taux de sensibilité ou "vrai positif"
table.BC[1,1]/sum(table.BC[1,]) # taux de spécificité ou de "vrai négatif"

precision_rate <- table.BC[2,2]/(table.BC[2,2] + table.BC[1,2])
precision_rate # taux de précision

confusion_afd <- table(y,pred.app$pred_class)
confusion_afd

# Taux d'erreur 
err_rate_afd <- 1 - sum(diag(confusion_afd)) / sum(confusion_afd)
err_rate_afd

# Taux de sensibilité
sens_rate_afd <- confusion_afd[2,2] / sum(confusion_afd[2,])
sens_rate_afd

# Taux de spécificité
spec_rate_afd <- confusion_afd[1,1] / sum(confusion_afd[1,])
spec_rate_afd

# Taux de précision
prec_rate_afd <- confusion_afd[2,2] / (confusion_afd[2,2] + confusion_afd[1,2])
prec_rate_afd

# Taux d'erreur en validation croisée
res.geoDA <- geoDA(X,y,validation="crossval")
res.geoDA$error_rate 

# Coefficients de corrélations avec les variables discriminantes
rap.cor <- res.desDA$values[1]/(1+res.desDA$values[1]) ##rapport de corrélation de la var discriminante
round(res.desDA$discor,3) #corrélation var discriminante et var X

FRatio(res.desDA$scores[,1],y) # F-Ratio composante 1
FRatio(res.desDA$scores[,2],y)# F-Ratio composante 2
FRatio(res.desDA$scores[,3],y)# F-Ratio composante 3
FRatio(res.desDA$scores[,4],y)# F-Ratio composante 4

corRatio(res.desDA$scores[,1],y) # corRatio composante 1
corRatio(res.desDA$scores[,2],y) # corRatio composante 2
corRatio(res.desDA$scores[,3],y) # corRatio composante 3
corRatio(res.desDA$scores[,4],y) # corRatio composante 4



# ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
# PLS DA AVEC FACTOMINER----
# ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
plDA_3 <- plsDA(X,y,autosel=FALSE,comps=3,cv="LOO")
plDA_3$R2 #rapport de corrélation
plDA_3$confusion
plDA_3$error_rate
plDA_3$Q2 #choix du nombre de composantes 
plot(plDA_3)


# Matrice de confusion
confusion_plsda <- plDA_3$confusion
confusion_plsda

# Taux d'erreur
err_rate_plsda <- 1 - sum(diag(confusion_plsda)) / sum(confusion_plsda)
err_rate_plsda

# Taux de sensibilité
sens_rate_plsda <- confusion_plsda[2,2] / sum(confusion_plsda[2,])
sens_rate_plsda

# Taux de spécificité
spec_rate_plsda <- confusion_plsda[1,1] / sum(confusion_plsda[1,])
spec_rate_plsda

# Taux de précision
prec_rate_plsda <- confusion_plsda[2,2] / (confusion_plsda[2,2] + confusion_plsda[1,2])
prec_rate_plsda

# Comparaison
plsda_res <- plsda(X, y, ncomp = 3)
perf.plsda <- perf(plsda_res, validation = "Mfold", folds = 5, 
                   progressBar = FALSE, auc = TRUE, nrepeat = 5)

# Affichage des courbes ROC et précision-rappel
plot(perf.plsda, ylim = c(0, 1), col = color.mixo(5:7), sd = TRUE, legend.position = "horizontal")

# Comparaison des matrices de confusion
confusion_lda
confusion_plsda


# Graphiques supplémentaires 

aa <- which(df_train[,1]=="aa")
ao <- which(df_train[,1]=="ao")
dcl <- which(df_train[,1]=="dcl")
iy <- which(df_train[,1]=="iy")
sh <- which(df_train[,1]=="sh")

plot(plDA_3$scores[,1],plDA_3$scores[,2], type='n')
points(plDA_3$scores[aa ,1],plDA_3$scores[aa ,2],pch=20,col="#FF0000")#rouge
points(plDA_3$scores[ao ,1],plDA_3$scores[ao ,2],pch=20,col="#FF9900")#orange
points(plDA_3$scores[dcl ,1],plDA_3$scores[dcl ,2],pch=20,col="#CC00CC")#violet
points(plDA_3$scores[iy ,1],plDA_3$scores[iy ,2],pch=20,col="#FFFF00")#jaune
points(plDA_3$scores[sh ,1],plDA_3$scores[sh ,2],pch=20,col="#3399FF")#bleu

# ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
# PLSDA AVEC MIXOMOICS----
# ''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
library(mixOmics)
library(RColorBrewer)
palette <- brewer.pal(5, "Set1") #palette de couleur pour les 5 composantes 

X_train <- as.matrix(df_train[,2:257])
Y_train <- as.factor(df_train[,1])
X_test <- as.matrix(df_test[,2:257])
Y_test <- df_test[, 1]

# On effectue PLSDA
plsda_1 <- plsda(X_train, Y_train, ncomp = 5)


# Utilisation de la fonction PERF 
perf <- perf(plsda_1, validation = "Mfold", folds = 5, 
                     progressBar = FALSE, auc = TRUE, nrepeat = 5) 

# Graphiques avec les différentes distances 
plot(perf, ylim=c(0,0.6),col = color.mixo(5:7), sd = TRUE, legend.position = "horizontal")

list.keepX <- c(1:5)  
tune.splsda_1 <- tune.splsda(X_train, Y_train, ncomp = 3, validation = 'Mfold', folds = 5,
                             progressBar = FALSE, dist = 'mahalanobis.dist', measure = "overall",
                             test.keepX = list.keepX, nrepeat = 10, cpus = 2)

# Nombre de composantes optimales 
ncomp <- tune.splsda_1$choice.ncomp$ncomp 
ncomp

error_1 <- tune.splsda_1$error.rate  # taux d'erreur par composantes
error_1

# PLSDA avec le nombre optimal de composantes
plsda_3 <- plsda(X_train, Y_train, ncomp = 3)

par(mfrow = c(1, 2))
# Représentation des individus PLSDA à 10 composantes
plotIndiv(plsda_1, group = Y_train, legend.title = "Phonème", ellipse = TRUE, legend = TRUE,
          col = palette,title = "PLS-DA sur 5 composantes")

# Représentation des individus PLSDA à 3 composantes
plotIndiv(plsda_3,, comp = c(1,3), group = Y_train, legend.title = "Phonème", ellipse = TRUE, legend = TRUE,
          title = "PLS-DA sur 3 composantes",col = palette)

pred_Y_test <- predict(plsda_3, newdata = X_test)$class
mean(pred_Y_test[["mahalanobis.dist"]] == df_test[,1]) # taux de classification
test.predict <- predict(plsda_3, X_test, dist = "max.dist")
Prediction <- test.predict$class$max.dist[, 3] 

# Matrice de confusion
confusion_PLSDA <- get.confusion_matrix(truth = Y_test, predicted = Prediction)
confusion_PLSDA

# Taux d'erreur
err_rate_plsda <- 1 - sum(diag(confusion_PLSDA)) / sum(confusion_PLSDA)
err_rate_plsda

# Taux de sensibilité
sens_rate_plsda <- confusion_PLSDA[2,2] / sum(confusion_PLSDA[2,])
sens_rate_plsda

# Taux de spécificité
spec_rate_plsda <- confusion_PLSDA[1,1] / sum(confusion_PLSDA[1,])
spec_rate_plsda

# Taux de précision
prec_rate_plsda <- confusion_PLSDA[2,2] / (confusion_PLSDA[2,2] + confusion_PLSDA[1,2])
prec_rate_plsda



