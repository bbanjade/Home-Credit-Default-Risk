
# Install Packages if not present
if(!require("tidyverse")) install.packages("tidyverse")
if(!require("ggplot2")) install.packages("ggplot2")
if(!require("dplyr")) install.packages("dplyr")
if(!require("e1071")) install.packages("e1071")
if(!require("caret")) install.packages("caret")
if(!require("zoo")) install.packages("zoo")
if(!require("reshape2")) install.packages("reshape2")
if(!require("plotly")) install.packages("plotly")
if(!require("data.table")) install.packages("data.table")
if(!require("gridExtra")) install.packages("gridExtra")
if(!require("MASS")) install.packages("MASS")
if(!require("PROC")) install.packages("PROC")
if(!require("fastAdaboost")) install.packages("fastAdaboost")
if(!require("cubist")) install.packages("cubist")
if(!require("C50")) install.packages("C50")


## Load libraries:
library(tidyverse)
library(dplyr)
library(ggplot2)
library(e1071)
library(caret)
library(zoo)
library(reshape2)
library(plotly)
library(data.table)
library(gridExtra)
library(MASS)
library(pROC)
library(fastAdaboost)
library(Cubist)
library(C50)


##################################################################################
###########Loading, cleaning, and preprocessing of data##############################
###################################################################################

data <- read.csv("BUAN6357_ProjectData_BANJADE.csv")
data <- data[, c(-1)]
str(data)
sum(is.na(data))

# filling up the na value with mean value of each column
data1 <- na.aggregate(data)
summary(data1)
sum(is.na(data1))

################################### skewness  #################################

numeric_list <- unlist(lapply(data1, is.numeric))
dt1_num <- setDT(data1)[,..numeric_list]
skewValues <- apply(dt1_num, na.rm = TRUE, 2, skewness)
skewValues

# BoxCox transformation
BoxCoxValues <- apply(dt1_num, 2, function(x) BoxCoxTrans(x, na.rm = TRUE))
x = list()
for (i in 1:ncol(dt1_num)){
  lambda <- BoxCoxValues[[i]][[1]]
  x[[i]] <- lambda
}

  # identification of the variables needed to be transformed
lambda = do.call(rbind, x)
lambda_df <- as.data.frame(cbind(colnames(dt1_num),lambda))
colnames(lambda_df)[1] <- "Column"
colnames(lambda_df)[2] <- "lambda"
knitr::kable(setDT(lambda_df)[!is.na(lambda)])

# transforming variables to reduce skewness
preProcValues <- preProcess(data1, method = "BoxCox")
preProcValues
dt1_tran <- predict(preProcValues, data1)

#Recreate numeric list with new dt1_tran
numeric_list <- unlist(lapply(dt1_tran, is.numeric))
dt1_num <- setDT(dt1_tran)[,..numeric_list]
dt1_num

# codes showing the effect of the transformation on dataset
col_trans <- lambda_df[!is.na(lambda)]$Column
i = 5
x <- list(
  title = as.character(col_trans[i])
)
p1 <- plot_ly(x = ~setDT(data)[,get(as.character(col_trans[i]))], type = "histogram", autobinx = FALSE) %>% layout(showlegend = FALSE) 
p2 <- plot_ly(x = ~setDT(dt1_tran)[,get(as.character(col_trans[i]))], type = "histogram", autobinx = FALSE) %>% layout(showlegend = FALSE)
subplot(p1,p2)


### creating doPlots and plotHist functions
doPlots <- function(data_in, fun, ii, ncol=3) {
  pp <- list()
  for (i in ii) {
    p <- fun(data_in=data_in, i=i)
    pp <- c(pp, list(p))
  }
  do.call("grid.arrange", c(pp, ncol=ncol))
}

plotHist <- function(data_in, i) {
  data <- data.frame(x=data_in[[i]])
  p <- ggplot(data=data, aes(x=x)) + geom_histogram(bins=100, fill="#0072B2", alpha = .9) + xlab(colnames(data_in)[i]) + theme_light() + 
    theme(axis.text.x = element_text(angle = 90, hjust =1))
  return (p)
}
# plots of the dataset before transformation
doPlots(as.data.frame(data1)[, (colnames(data1) %in% 
                                  as.character(col_trans))], plotHist, ii = 1:length(col_trans))
# plots of the dataset after transformation
doPlots(as.data.frame(dt1_tran)[, (colnames(dt1_tran) %in% 
                                     as.character(col_trans))], plotHist, ii = 1:length(col_trans))


################ correlation coefficient#############################

data11 <- dt1_num[, c(5:13,15:18)] ## only numeric values are considered
names(data11)
cor(data11, method = c("pearson", "kendall", "spearman"))
  
#compute the correlation matrix
cormat <- round(cor(data11),2)
melted_cormat <- melt(cormat)
head(melted_cormat)
ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile()

# correlation matrix has redundant information. Function are used below to set half of it to NA.
# Get lower triangle of the correlation matrix
get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}
# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

upper_tri <- get_upper_tri(cormat)
upper_tri

# Melt the correlation data and drop the rows with NA values
# Melt the correlation matrix
melted_cormat <- melt(upper_tri, na.rm = TRUE)
ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()

## Reorder the correlation matrix
## Helper function to reorder the correlation matrix:
reorder_cormat <- function(cormat){
  # Use correlation between variables as distance
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
}

cormat <- reorder_cormat(cormat)
upper_tri <- get_upper_tri(cormat)
# Melt the correlation matrix
melted_cormat <- melt(upper_tri, na.rm = TRUE)
# Create a ggheatmap
ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()
# Print the heatmap
print(ggheatmap)

# boxplots to know the outliers
boxplot(dt1_num$REGION_RATING_CLIENT)
boxplot(dt1_num$REGION_RATING_CLIENT_W_CITY)
boxplot(dt1_num$AMT_ANNUITY)
boxplot(dt1_num$AMT_GOODS_PRICE)
boxplot(dt1_num$AMT_CREDIT)

# removal of highly correlated variables
data12 <- dt1_num[, -c(5,8,9,12)] 
names(data12)

############################## Descriptive Statistics #################################
###########################################################################

data12 %>%group_by(TARGET)%>%summarize(n()) # 0 are 18302 and 1 are 1697

ggplot(data = melt(data12), mapping = aes(x = value)) + 
  geom_histogram(bins = 30) + facet_wrap(~variable, scales = 'free_x')


############################################################################################
##############################     MODELING PART    ########################################
###########################################################################################

# preparation for model
data12_prep <- mutate(data12, TARGET = ifelse(TARGET == 0,"Yes","No"))
data12_prep$TARGET <- as.factor(data12_prep$TARGET)

training.samples <- createDataPartition(data12_prep$TARGET, p = .8)[[1]]
train.data <- data12_prep[ training.samples, ]
test.data  <- data12_prep[-training.samples, ]
View(train.data)

traincntrl <- trainControl(method = 'repeatedcv',
                           number = 5,
                           repeats = 2,
                           classProbs = TRUE, 
                           sampling = "down",
                           summaryFunction = twoClassSummary)

## Logisitc Regression
## step logistic regression model
model <- glm(TARGET~., data = train.data, family = binomial) %>% stepAIC(trace = FALSE)
summary(model) 
#Make Predictions
probabilities <- model %>% predict(test.data, type = "response")
predicted.classes <- ifelse(probabilities > 0.5, "No", "Yes")
#Model accuracy
mean(predicted.classes==test.data$TARGET)

## Full logistic regression model
full.model <- glm(TARGET ~., data = train.data, family = binomial)
summary(full.model)
# Make predictions
probabilities <- predict(step.model, test.data, type = "response")
predicted.classes <- ifelse(probabilities > 0.5, "No", "Yes")
# Prediction accuracy
observed.classes <- test.data$TARGET
mean(predicted.classes == observed.classes)


# logistic regression for project
logisticReg <- train(TARGET ~.,
                     data = train.data,
                     family = binomial,
                     trControl = traincntrl)

logisticReg
plot(logisticReg, scales = list(x=list(log(2))))

# 2.2) Running a KNN model
knnFit <- train(TARGET ~.,
                data = train.data,
                method = "knn",
                preProc = c("center", "scale"),
                tuneGrid = data.frame(.k = 3:6),
                trControl = traincntrl)
knnFit
plot(knnFit, scales = list(x=list(log=2)))

# 2.3) Running a SVM model
svmFit <- train(TARGET ~.,
                data = train.data,
                method = 'svmRadial',
                preProc = c('center','scale'),
                tuneLength = 7,
                trControl = traincntrl)

svmFitRadial <- svmFit
svmFitLinear <- train(TARGET ~.,
                      data = train.data,
                      method = 'svmLinear',
                      preProc = c('center','scale'),
                      #metric = "ROC",
                      tuneLength = 7,
                      trControl = traincntrl)

#### Comparing all the SVM models
resamp <- resamples(list(SVM_Radial = svmFit, SVM_Linear = svmFitLinear))
summary(resamp)

## Comparing SVM, Logistic Regression and KNN

resamp <- resamples(list(SVM = svmFitLinear, Logistic = logisticReg, KNN = knnFit))
summary(resamp)
resamp


#  Measuring Performance in Classification Models
# Class Predictions

test.data$svmFitLinearclass <- predict(svmFitLinear, test.data)
test.data$svmFitLinearprobs <- predict(svmFitLinear, newdata = test.data , type = "prob")

test.data$logclass <- predict(logisticReg, test.data)
test.data$logprobs <- predict(logisticReg, newdata = test.data , type = "prob")

test.data$knnFitclass <- predict(knnFit, test.data)
test.data$knnFitprobs <- predict(knnFit, newdata = test.data , type = "prob")

## Evaluating Predicted Classes
confusionMatrix(data = test.data$svmFitLinearclass,
                reference = test.data$TARGET,
                positive = "Yes")

## Following is the confusion matrix of the Logistic Regression
confusionMatrix(data = test.data$logclass,
                reference = test.data$TARGET,
                positive = "Yes")
######## for knn
confusionMatrix(data = test.data$knnFitclass, 
                reference = test.data$TARGET, 
                positive = "Yes")


## Evaluating Class Probabilities
rocCurve <- roc(response = test.data$TARGET, predictor = test.data$svmFitLinearprobs[,1], levels = rev(levels(test.data$TARGET)))
plot(rocCurve, legacy.axes = TRUE)
auc(rocCurve)
## Area under the curve: 0.7246..

## Logistic Regression
rocCurve <- roc(response = test.data$TARGET, predictor = test.data$logprobs[,1], levels = rev(levels(test.data$TARGET)))
plot(rocCurve, legacy.axes = TRUE)
auc(rocCurve)
## Area under the cruve: 0.7173.

## knn
rocCurve <- roc(response = test.data$TARGET, predictor = test.data$knnFitprobs[,1], levels = rev(levels(test.data$TARGET)))
plot(rocCurve, legacy.axes = TRUE)
auc(rocCurve)
## Area under the cruve is 0.6627

########################## Classification Trees  ##################
########################################################################

dtFitCART <- train(x= setDT(train.data)[,-c('TARGET')], y= train.data$TARGET, method = 'rpart',
                   preProc = c('center','scale'), tuneLength = 7, metric = "ROC", trControl = traincntrl)
dtFitBagged <- train(x= setDT(train.data)[,-c('TARGET')], y= train.data$TARGET, method = 'treebag', preProc = c('center','scale'),
                     tuneLength = 7, metric = "ROC", trControl = traincntrl)
dtFitrf <- train(x= setDT(train.data)[,-c('TARGET')], y= train.data$TARGET, method = 'rf', preProc = c('center','scale'), tuneLength = 7,
                 metric = "ROC", trControl = traincntrl)


dtFitAdaboost <- train(x= setDT(train.data)[,-c('TARGET')],
                       y= train.data$TARGET,
                       method = 'adaboost',
                       preProc = c('center','scale'),
                       tuneLength = 7,
                       metric = "ROC",
                       trControl = traincntrl)

dtFitC5.0 <- train(x= setDT(train.data)[,-c('TARGET')],
                   y= train.data$TARGET,
                   method = 'C5.0',
                   preProc = c('center','scale'),
                   tuneLength = 7,
                   metric = "ROC",
                   trControl = traincntrl)
dtFitC5.0
alltreemodels <- resamples(list(CART = dtFitCART,  
                                Bagged = dtFitBagged,
                                RF = dtFitrf, AdaBoost = dtFitAdaboost, 
                                C5.0 = dtFitC5.0))
summary(alltreemodels)

plot(dtFitC5.0)

test.data$C5.0class <- predict(dtFitC5.0, test.data)
test.data$C5.0probs <- predict(dtFitC5.0, newdata = test.data , type = "prob")

confusionMatrix(data = test.data$C5.0class,
                reference = test.data$TARGET,
                positive = "Yes")

rocCurve <- roc(response = test.data$TARGET, predictor = test.data$C5.0probs[,1], levels = rev(levels(test.data$TARGET)))
plot(rocCurve, legacy.axes = TRUE)
auc(rocCurve)
