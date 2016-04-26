library(dplyr)
library(readr)
library(ggplot2)
library(caret)
library(pROC)
library(ROCR)
library(parallel)
library(randomForest)
library(lubridate)

###############################################################################
# Function Name: read_yactives (Helper Function)
# Inputs: file_name: A string representing the file name of the yactives file.
#         
# Returns a data frame containing entire yactives data set.
###############################################################################
read_yactives <- function(file_name){
  yactives <- read_csv("yactives_training_set_June10_A.csv")
  
  return(yactives)
}

###############################################################################
# Function Name: get_classes (Helper Function)
# Inputs: yactives: A data frame containing features and outcome classes.
#                   Outcome classes are assumed to be in the form P1, P2, 
#                   P3, ... .  The drug ID column is also removed from the data.
#         
# Returns a data frame containing classes only.
###############################################################################
get_classes <- function(yactives){
  # The following transformation is needed to make the class data 
  # compatible with the training and prediction functions.
  # Convert 1's to true and 0's to false in the class columns
  # Then, convert the true and false values to factors.  In order
  # to calculate the probabilities of each class value, the class
  # value has to be a legal variable name (1 and 0 are not legal 
  # variable names).  The classes also need to be factors in 
  # order to be compatible with train and predict.  Note that the 
  # class columns are P1 through P10
  classes <- yactives %>% 
    select(starts_with("P")) %>%
    mutate_each(funs(ifelse(. == 1, "true", "false"))) %>%
    mutate_each(funs(as.factor(.)))
  
  return(classes)
}

###############################################################################
# Function Name: get_features (Helper Function)
# Inputs: yactives: A data frame containing features and outcome classes.
#                   Outcome classes are assumed to be in the form P1, P2, 
#                   P3, ... .  The drug ID column is also removed from the data.
#         
# Returns a data frame containing features only.
###############################################################################
get_features <- function(yactives){
  features <- select(yactives, -starts_with("P"), -ID)
  return(features)
}

###############################################################################
# Function Name: create_train_test_sets (Helper Function)
# Inputs: yactives: The entire yactives data set.
#         p_train: The proportion of training data out of the entire 
#                  yactives training data set.
#         
# Returns a list containing training data, training class outcomes,
# test data, test class outcomes, and phenotype class names 
# (like P1, P2, ...).  Elements of the list can be accessed with list 
# tag values train_class, train_set, test_class, test_set, and class_names.
# For example, assuming the return value was named "test_train_data", access 
# the training data set as follows: test_train_data$test_set or
# test_train_data[["test_set"]].
###############################################################################
create_train_test_sets <- function(yactives, p_train){
  # create a new class which is the logical OR of P7 and P8.  P7 and P8
  # represent zebra fish death at 24 and 48 hours.  Combining these 2 
  # phenotypes creates a single death predictor.
  yactives <- yactives %>% mutate(P11 = P7 | P8)
  
  # get the indices for the training partition
  training_part <- createDataPartition(y = yactives$P1, p = p_train)
  
  # create training and test data frames
  train_set <- slice(yactives, training_part$Resample1)
  test_set <- slice(yactives, -training_part$Resample1)
  
  # separate the train phenotype classes from the features
  train_class <- get_classes(train_set)
  train_set <- get_features(train_set)
  
  # separate the test phenotype classes from the features
  test_class <- get_classes(test_set)
  test_set <- get_features(test_set)
  
  # if the small dataset flag is true
  if (small_yactives_test_data_enabled){
    # reduce the data size
    test_class <- select_(test_class, .dots = small_test_class)
    train_class <- select_(train_class, .dots = small_test_class)
  }
  
  # create a vector of class names
  class_names <- colnames(test_class)
  
  # return a list of class and feature data frames
  return_list <- 
    list(train_class = train_class, 
         train_set = train_set, 
         test_class = test_class, 
         test_set = test_set, 
         class_names = class_names)
  return(return_list)
}

###############################################################################
# Function Name: process_knn (process KNN model)
# Description: Performs all processing related to generating the KNN model 
#             for a specific training feature set and outcome class.
# Inputs: train_set: A training set data frame excluding outcome classes.
#         train_class: Class outcome single column data frame for the train set.
#         test_set: A test set data frame excluding outcome classes.
#         test_class: Class outcome single column data frame for the test set.
#         num_cross_validations: The number of cross validations to perform.
#         train_percentage: The proportion of test_set and test_class to use
#                           for training.
#         
# Returns a list containing fit, predict, roc objects.  Elements of the list 
# can be accessed with list tag values "fit", prdict", and "roc"
###############################################################################
process_knn <- function(train_set, train_class, 
                        test_set, test_class,
                        num_cross_validations = 10,
                        train_percentage = 0.8,
                        max_k = 0){
  # create a training control object which includes data to enable ROC summary
  control <- trainControl(method = 'cv', 
                          number = num_cross_validations, 
                          p = train_percentage, 
                          classProbs = TRUE, summaryFunction = twoClassSummary)
  
  # if max_k was not overridden by the caller
  if (max_k == 0){
    # calculate the max K
    # wew to do: review this
    max_k <- length(train_class) / 2
    # if max_k is even
    if (max_k %% 2 == 0){
      # make max_k odd
      max_k <- max_k - 1
    }
  }
  
  knn_fit <- train(train_set, train_class, method = "knn", 
                   trControl = control, tuneLength = 1, 
                   tuneGrid = data.frame(k = seq(3, max_k, 2)))
  
  # fit the test set to the KNN model and create a prediction
  knn_predict <- predict(knn_fit, test_set, type = "prob")
  
  #Get the confusion matrix to see accuracy value and other parameter values
  # wew to do: Fix this
  # confMatrix <- confusionMatrix(knn_predict, test_set$P1)
  # confMatrix
  
  # plot the ROC curve
  knn_roc <- roc(test_class, knn_predict[,"true"], 
                 levels = rev(test_class))
  plot(knn_roc)
  
  ret_list <- list(fit = knn_fit, predict = knn_predict, roc = knn_roc)
  return(ret_list)
}

###############################################################################
# Function Name: process_rf (process random forest model)
# Description: Performs all processing related to generating the random
#              forest model for a specific training feature set and outcome 
#              class.
# Inputs: train_set: A training set data frame excluding outcome classes.
#         train_class: Class outcome single column data frame for the train set.
#         test_set: A test set data frame excluding outcome classes.
#         test_class: Class outcome single column data frame for the test set.
#         num_cross_validations: The number of cross validations to perform.
#         train_percentage: The proportion of test_set and test_class to use
#                           for training.
#         
# Returns a list containing fit, predict, roc objects.  Elements of the list 
# can be accessed with list tag values "fit", prdict", and "roc".
###############################################################################
process_rf <- function(train_set, train_class, 
                       test_set, test_class){
  # fit the test set to the random forest model
  rf_fit <- randomForest(train_set, train_class, importance = TRUE, proximity = TRUE)
  
  # create the prediction object
  rf_predict <- prediction(rf_fit$votes[,2], train_class)
  
  #AUC score
  auc_score <- performance(prediction.obj = rf_predict, measure = 'auc')@y.values
  
  # package the fit, prediction, and score into a list and return it to caller
  ret_list <- list(fit = rf_fit, predict = rf_predict, auc_score = auc_score)
  return(ret_list)
}

###############################################################################
# Function Name: process_models
# Description: Top level function which coordinates collecting model data
#              for all model types lik KNN, random forest, etc.
# Inputs: train_test_list: A list containing sublists of training factors, 
#                          training classes, test factors, test classes,
#                          and class names (P1, P2, P3, ...)
# Returns a list containing sublists of model type, phenotype, and model summary 
# For example, if the return value is assigned to the variable 
# "model_results_list",  access the P3 KNN ROC object as follows: 
# model_results_list$knn$P3$roc
###############################################################################
process_models <- function(train_test_list){
  # extract the individual train and test frames out of the master list
  class_names <- train_test_list[["class_names"]]
  train_set <- train_test_list[["train_set"]]
  train_class <- train_test_list[["train_class"]]
  test_set <- train_test_list[["test_set"]]
  test_class <- train_test_list[["test_class"]]
  
  # create empty lists for each model type
  knn_list <- list()
  rf_list <- list()
  
  # for each phenotype class in the class list
  for(phen_class in class_names){
    # extract the training class vector from the training class frame
    train_class_n <- select_(train_class, .dots = phen_class)
    train_class_n <- train_class_n[[1]]
    
    # extract the test class vector from the training class frame
    test_class_n <- select_(test_class, .dots = phen_class)
    test_class_n <- test_class_n[[1]]  
    
    # if knn model processing is enabled
    if (knn_enabled){
      cat("Running KNN Model for class", phen_class, "\n")
      
      # run the KNN model
      start_time <- Sys.time()
      knn_results_list <- 
        process_knn(train_set, train_class_n, test_set, test_class_n)
      duration <- difftime(Sys.time(), start_time, units = "mins")
      cat("Run Time:", as.numeric(duration), "minutes\n")
      
      # save the model results in the KNN model list
      knn_list[[phen_class]] <- knn_results_list
    }
    
    # if random forest processing is enabled
    if (rf_enabled){
      cat("Running Random Forest Model for class", phen_class, "\n")
      
      # run the random forest model
      start_time <- Sys.time()
      rf_results_list <- 
        process_rf(train_set, train_class_n, test_set, test_class_n)
      duration <- difftime(Sys.time(), start_time, units = "mins")
      cat("Run Time:", as.numeric(duration), "minutes\n")
      
      # save the model results in the random forest model list
      rf_list[[phen_class]] <- rf_results_list  
    }
  }
  
  # concatenate the model lists together and return to the caller
  ret_model_results_list <- list()
  ret_model_results_list$knn <- knn_list
  ret_model_results_list$rf <- rf_list
  return(ret_model_results_list)
}

###############################################################################
# main program starts here
###############################################################################
# use the following variables to enable and disable specific models.
knn_enabled = FALSE
rf_enabled = TRUE

# use the following variables to reduce the data size by selecting a single 
# phenotype class.
small_yactives_test_data_enabled <- FALSE
small_test_class <- "P11"

# read the training data file
yactives <- read_yactives("yactives_training_set_June10_A.csv")

# get the training and test sets
train_prop <- 0.8
train_test_list <- create_train_test_sets(yactives, train_prop)

# Train the models
# Notes on how to use the model_results_list.  The model results 
# list is a list of lists of lists. The top level is a list of models 
# like knn and rf.  The second level is a list phenotypes like P1, P2, 
# P3, ... .  The third level is a list of model summary parameters like 
# training object, fit object, roc object.  To access model parameters, 
# one can use the '$' notation.  For example, access the KNN P3 ROC 
# model as follows: model_results_list$knn$P3$roc
model_results_list <- process_models(train_test_list)

# print the ROC AUC results for KNN
auc_results <- lapply(model_results_list$rf, function(i){
  i$auc_score
})
auc_results