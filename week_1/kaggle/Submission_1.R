
# Download 1_random_forest_r_submission.csv from the output below
# and submit it through https://www.kaggle.com/c/titanic-gettingStarted/submissions/attach
# to enter this getting started competition!

library(ggplot2)
library(randomForest)
library(mgcv)
library(MASS)
library(gbm)
library(klaR)
library(e1071)
library(FNN)

set.seed(1)
setwd("/Users/jamesdiao/Documents/Kohane_Lab/kaggle")
train_orig <- read.csv("./train.csv", stringsAsFactors=FALSE)
valid_flag <- runif(nrow(train_orig)) > 0.7
train <- train_orig[valid_flag,]
test <- train_orig[!valid_flag,]
train_y <- train$Survived
test_y <- test$Survived
test_orig  <- read.csv("./test.csv",  stringsAsFactors=FALSE)

extractFeatures <- function(data) {
  features <- c("Pclass",
                "Age",
                "Sex",
                "Parch",
                "SibSp",
                "Fare",
                "Embarked")
  fea <- data[,features]
  fea$Age[is.na(fea$Age)] <- -1
  fea$Fare[is.na(fea$Fare)] <- median(fea$Fare, na.rm=TRUE)
  fea$Embarked[fea$Embarked==""] = "S"
  fea$Sex      <- as.factor(fea$Sex)
  fea$Embarked <- as.numeric(fea$Embarked == "C")
  fea$Sex      <- as.numeric(fea$Sex)-1
  return(fea)
}

train <- extractFeatures(train)
test <- extractFeatures(test)
train <- cbind(train_y, train)

str(train)
#gam_init <- gam(train[,1] ~ s(Age)+s(Fare), data = train)
#gam_init <- gam(train$Survived ~ extractFeatures(train))
#gam_y <- predict(gam_init, test)
#mean(gam_y > 0.5)

lm_init <- lm(train_y ~ ., data = train)
lm_y <- predict.lm(lm_init,test)
mean((lm_y > 0.5) == test_y)

knn_init <- knn(train[,-1], test, cl = train_y, k=50)
x <- as.numeric(knn_init[1:618])-1
mean(x == test_y)

logit_init <- glm(train_y ~ ., data = train, family = binomial("logit"))
logit_y <- predict.glm(logit_init, test)
mean((logit_y > 0.5) == test_y)

svm_init <- svm(train_y ~ ., data = train, kernel = "linear")
svm_y <- predict(svm_init,test)
mean((svm_y > 0.5) == test_y)

nb_init <- NaiveBayes(as.factor(train_y)~ ., data = train)
nb_y <- predict(nb_init,test)
if (is.factor(nb_y$class)) nb_y <- as.numeric(nb_y$class)-1
mean(nb_y == test_y)

gbm_init <- gbm.fit(train[,-1], train_y, distribution = "bernoulli", n.trees = 2000L, shrinkage = 0.001, interaction.depth=3)
plot(gbm_init)
summary(gbm_init, plotit = T)
gbm_y <- predict(gbm_init,extractFeatures(test), n.trees = 2000L) 
mean((gbm_y > 0.5) == test_y)

rf <- randomForest(train[,-1], as.factor(train_y), ntree=1, importance=TRUE)
plot(rf)
rf_y <- as.numeric(predict(rf, test))-1
mean(rf_y == test_y)
imp <- importance(rf, type=1)
imp
plot(imp)

#test$lm_pred <- lm_y
#train$lm_pred <- predict(lm_init,train, n.trees = 2000L)
test$gbm_pred <- gbm_y
train$gbm_pred <- predict(gbm_init,train, n.trees = 2000L)
temp <- predict(nb_init,test)
test$nb_pred <- temp$posterior[,2]
temp <- predict(nb_init,train)
train$nb_pred <- temp$posterior[,2]

surv_final <- svm_y
submission <- data.frame(PassengerId = test$PassengerId)
submission$Survived <- surv_final
mean(submission$Survived == test_y)

write.csv(submission, file = "submission_1.csv", row.names=FALSE)










train$Survived <- factor(train$Survived, levels=c(1,0))
levels(train$Survived) <- c("Survived", "Died")
train$Pclass <- as.factor(train$Pclass)
levels(train$Pclass) <- c("1st Class", "2nd Class", "3rd Class")
train$Gender <- factor(train$Sex, levels=c("female", "male"))
levels(train$Gender) <- c("Female", "Male")
#png("1_survival_by_class.png", width=800, height=600)
mosaicplot(train$Pclass ~ train$Survived, main="Passenger Survival by Class",
           color=c("#8dd3c7", "#fb8072"), shade=FALSE,  xlab="", ylab="",
           off=c(0), cex.axis=1.4)
#dev.off()
#png("2_survival_by_gender.png", width=800, height=600)
mosaicplot(train$Gender ~ train$Survived, main="Passenger Survival by Gender",
           color=c("#8dd3c7", "#fb8072"), shade=FALSE,  xlab="", ylab="",
           off=c(0), cex.axis=1.4)
#dev.off()