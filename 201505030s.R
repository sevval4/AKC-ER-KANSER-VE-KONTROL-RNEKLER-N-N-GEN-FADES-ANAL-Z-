library(BiocManager)
library(GEOquery)
okunan1 <- getGEO("GSE18842", GSEMatrix = TRUE)
eset1 <- okunan1[[1]]
kopya <- exprs(eset1)
dim(eset1)
head(eset1)


BiocManager::install("hgu133plus2.db", force = TRUE)
annotation(eset1) <- "hgu133plus2.db"
colnames(pData(eset1))
durum <- factor(pData(eset1)$source_name_ch1)

library(genefilter)
# Veri filtreleme
filtrelenmis <- varFilter(eset1, var.cutoff = 0.9)
sonveri <- data.frame(t(exprs(filtrelenmis)))
sonveri$durum <- durum
dim(sonveri)

library(caret)
#  eğitim ve test setlerine ayırma
set.seed(111)
train_indices <- createDataPartition(durum, p = 0.7, list = FALSE)
train_data <- sonveri[train_indices, ]
test_data <- sonveri[-train_indices, ]

# RANDOM FOREST
model_rf <- train(durum ~ ., data = train_data, method = "rf")
predictions_rf <- predict(model_rf, newdata = test_data)
conf_matrix_rf <- confusionMatrix(predictions_rf, test_data$durum)
accuracy_rf <- conf_matrix_rf$overall["Accuracy"]
print(paste("Random Forest Doğruluk DEĞERİ:", accuracy_rf))
feature_importance_rf <- varImp(model_rf)
print(feature_importance_rf)

#SVM
svm_model <- train(durum ~ ., data = train_data, method = "svmLinear")
predictions_svm <- predict(svm_model, newdata = test_data)
conf_matrix_svm <- confusionMatrix(predictions_svm, test_data$durum)
accuracy_svm <- conf_matrix_svm$overall["Accuracy"]
print(paste("SVM Doğruluk DEĞERİ:", accuracy_svm))
feature_importance_svm <- varImp(svm_model)
print(feature_importance_svm)
print(conf_matrix_svm)


library(rpart)
# KARAR AĞACI
tree_model <- rpart(durum ~ ., data = train_data, method ="class")
predictions_tree <- predict(tree_model, newdata = test_data, type = "class")
conf_matrix_tree <- confusionMatrix(predictions_tree, test_data$durum)
accuracy_tree <- conf_matrix_tree$overall["Accuracy"]
print(paste("Karar Ağacı Doğruluk:", accuracy_tree))
importance_tree <- tree_model$variable.importance
top_features <- importance_tree[order(-importance_tree, decreasing = TRUE)][1:20]
print(top_features)
print(conf_matrix_tree)