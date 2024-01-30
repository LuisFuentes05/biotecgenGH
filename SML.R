library(tidyverse)
library(ggplot2)
library(plotly)
library(FactoMineR)

normado <- read.table('/Users/luisfuentes/Desktop/biotecgen/newData/normandeGenotypedCodeIlu.txt', header=TRUE, sep='\t', na.strings=c(""))
normado[is.na(normado)] <- 20

rsByFilt <- read.table('/Users/luisfuentes/Desktop/biotecgen/selectedRSTraitDB.txt', header=TRUE, sep='\t')

rsByFiltT <- t(rsByFilt)
colnames(rsByFiltT) <- rsByFilt$rs_ID
rsByFiltT <- rsByFiltT[-1,]
pca_result <- prcomp(normado[, -1], center=TRUE, scale.=TRUE)

normadoFil <- merge(normado, rsByFilt, by.x='dbSNP_RS_ID', by.y='rs_ID', all.x=FALSE)

normadoFilC <- normadoFil
normadoFilC <- normadoFilC[, !(names(normadoFilC) %in% c('dbSNP_RS_ID', 'Affy_SNP_ID'))]
rownames(normadoFilC) <- normadoFilC$AX_Id
normadoFilC <- normadoFilC[, -which(names(normadoFilC) %in% c('AX_Id'))]

# Transposing the dataframe
normadoFilCT <- t(normadoFilC)
colnames(normadoFilCT) <- NULL
normadoFilCT <- data.frame(apply(normadoFilCT, 2, as.integer))

# Calculating threshold for dropping columns
threshold <- ceiling(nrow(normadoFilCT) * 0.9) + 1

normadoFilCT <- normadoFilCT[, -which(colSums(normadoFilCT == 20) >= nrow(normadoFilCT) - threshold)]

# Describing the data and filtering based on standard deviation
fil <- as.data.frame(summary(normadoFilCT))
fil <- subset(fil, std > 0 & std < 1)[, -c('count', 'mean', 'min', '25%', '50%', '75%', 'max')]

newFilAgr <- normadoFilCT[, fil$rs]

newFilAgr <- t(newFilAgr)
colnames(newFilAgr) <- newFilAgr[1, ]
newFilAgr <- newFilAgr[-1, ]

# Calculating correlation
corr <- cor(newFilAgr, method = "pearson")

# Counting negative and positive correlations
corrRSN <- colSums(corr < -0.7)
corrRSP <- colSums(corr > 0.7)

# Filtering columns based on conditions
RS2fil <- names(corr)[corrRSP > 1 & corrRSN > 0]
secFilAgr <- newFilAgr[, !(colnames(newFilAgr) %in% RS2fil)]

# Plotting the heatmap
heatmap(corr, main="Matriz de Correlación entre variables", col=cm.colors(256), scale="none", symm=TRUE)


#### SUPERVISED MACHINE LEARNING
normadoScored2$Score <- ifelse(normadoScored2$Score == 20, -1, normadoScored2$Score)

Y <- normadoScored2$Score
X <- normadoScored2[, 2:437]

# Splitting the data into training and testing sets
set.seed(41)
index <- sample(1:nrow(X), 0.8 * nrow(X))
x_train <- X[index, ]
x_test <- X[-index, ]
y_train <- Y[index]
y_test <- Y[-index]

# Linear Regression
library(Metrics)

clf <- lm(y_train ~ ., data = cbind(y_train, x_train))
y_pred <- predict(clf, newdata = x_test)
mse <- mse(y_test, y_pred)
mae <- mae(y_test, y_pred)

cat('Mean Absolute Error (MAE):', mae, '\n')
cat('Mean Squared Error (MSE):', mse, '\n')

# Logistic Regression
library(nnet)

model <- multinom(as.factor(y_train) ~ ., data = cbind(y_train, x_train))
predictions <- predict(model, newdata = x_test, type = "response")
predicted_classes <- max.col(predictions) - 1
accuracy <- mean(predicted_classes == y_test)

cat('Accuracy:', accuracy, '\n')

# CatBoost Classifier
library(Catboost)

modelCat <- catboost.train(
  data = catboost.load_pool(data = x_train, label = y_train),
  params = list(
    depth = c(7, 8, 10),
    learning_rate = c(0.01, 0.05, 0.1),
    iterations = c(100, 200)
  ),
  nfold = 2,
  loss_function = 'Logloss'
)

best_params <- catboost.get_params(modelCat)
preds2 <- predict(modelCat, newdata = x_test)
accuracy2 <- mean(round(preds2) == y_test)

cat('Best Parameters:', best_params, '\n')
cat('Accuracy:', accuracy2, '\n')

# CatBoost Regressor
modelCatReg <- catboost.train(
  data = catboost.load_pool(data = x_train, label = y_train),
  params = best_params,
  nfold = 2,
  loss_function = 'RMSE'
)

best_estimator <- catboost.get_params(modelCatReg)
preds <- predict(modelCatReg, newdata = x_test)
mae <- mae(y_test, preds)
mse <- mse(y_test, preds)

cat('Best Estimator:', best_estimator, '\n')
cat('Mean Absolute Error (MAE):', mae, '\n')
cat('Mean Squared Error (MSE):', mse, '\n')
