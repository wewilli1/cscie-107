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
# Function Name: create_train_test_indep_sets (Helper Function)
# Inputs: yactives: The entire yactives data set.
#         train_prop: The proportion of training data out of the entire 
#                     yactives training data set.
#         yactives_indep: An independent yacives drug data set (cancer drugs)
#         
# Returns a list containing training data, training class outcomes,
# test data, test class outcomes, phenotype class names 
# (like P1, P2, ...), and drug ID's.  Elements of the list can be accessed 
# with list tag values.  For example, assuming the return value was named 
# "test_train_data", access the training data set as follows: 
# test_train_data$test_set or test_train_data[["test_set"]].
###############################################################################
create_train_test_indep_sets <- 
  function(yactives_train, train_prop, yactives_indep){
  # P2, P7, and P8 are "death" classes.  Create combinations of death 
  # classes by logically ORing the permutations of P2, P7, and P8
  yactives_train <- yactives_train %>% 
    mutate(P11 = P2 | P7 | P8, P12 = P2 | P7, P13 = P2 | P8, P14 = P7 | P8)
  
  # save the entire training set and class for the random forest (RF) model
  # because we don't have to split into train / test in the RF model
  train_set_rf <- get_features(yactives_train)
  train_class_rf <- get_classes(yactives_train)
  train_set_drug_id_rf <- select(yactives_train, ID)
  
  # For the rest of the models, create train and test by splitting into 
  # train / test sets.  Start by getting the indices for the training partition
  training_part <- createDataPartition(y = yactives_train$P1, p = train_prop)
  
  # create training and test data frames
  train_set <- slice(yactives_train, training_part$Resample1)
  test_set <- slice(yactives_train, -training_part$Resample1)
  
  # Save the drug ID's for the train and test sets
  train_set_drug_id <- select(train_set, ID)
  test_set_drug_id <- select(test_set, ID)
  
  # separate the train phenotype classes from the features
  train_class <- get_classes(train_set)
  train_set <- get_features(train_set)
  
  # separate the test phenotype classes from the features
  test_class <- get_classes(test_set)
  test_set <- get_features(test_set)
  
  # if the small dataset flag is true
  if (small_yactives_test_data_enabled){
    # reduce the data size by selecting a single phenotype class
    test_class <- select_(test_class, .dots = small_test_class)
    train_class <- select_(train_class, .dots = small_test_class)
  }
  
  # create a vector of class names
  class_names <- colnames(test_class)
  
  # Remove the duplicate "ID" column from the yactives_indep data
  yactives_indep <- yactives_indep[,-2]
  
  # save the yactives_indep drug ID's
  indep_set_drug_id <- select(yactives_indep, ID)
  
  # remove the drug ID's from yactives_indep
  yactives_indep <- select(yactives_indep, -ID)
  
  # return a list of class and feature data frames
  return_list <- 
    list(# random forest training set
         train_set_rf = train_set_rf,
         train_set_drug_id_rf = train_set_drug_id_rf,
         train_class_rf = train_class_rf,
         # split training set
         train_set = train_set,
         train_set_drug_id = train_set_drug_id,
         train_class = train_class,
         # split test set
         test_set = test_set,
         test_set_drug_id = test_set_drug_id,
         test_class = test_class, 
         # independent set
         indep_set = yactives_indep,
         indep_set_drug_id = indep_set_drug_id,
         # prediction class names
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
#         indep_set: The independent cancer drug data set
#         
# Returns a list containing fit, predict, roc objects.  Elements of the list 
# can be accessed with list tag values "fit", prdict", and "roc".
###############################################################################
process_rf <- function(train_set, train_class, indep_set){
  # fit the test set to the random forest model
  rf_fit <- randomForest(train_set, train_class, importance = TRUE, proximity = TRUE)
  
  # create the prediction object
  oob_vote_col <- 2
  train_predict <- prediction(rf_fit$votes[,oob_vote_col], train_class)
  
  # use the model to make a prediction on the independent data set
  indep_predict <- predict(rf_fit, type = "prob", newdata = indep_set)
  
  #AUC score
  # wew to do: Only call performance on the test set.
  auc_score <- performance(prediction.obj = train_predict, measure = 'auc')@y.values
  
  # package the fit, prediction, and score into a list and return it to caller
  ret_list <- list(fit = rf_fit, train_predict = train_predict, 
                   indep_predict = indep_predict, auc_score = auc_score)
  
  # return the random forest model list to the caller
  return(ret_list)
}

###############################################################################
# Function Name: process_models
# Description: Top level function which coordinates collecting model data
#              for all model types like KNN, random forest, etc.  Each model
#              can be enabled or disabled independently using a global 
#              flag located in main.
# Inputs: train_test_list: A list containing sublists of training factors, 
#                          training classes, test factors, test classes,
#                          and class names (P1, P2, P3, ...)
# Returns a list containing sublists of model type, phenotype, and model summary 
# For example, if the return value is assigned to the variable 
# "model_results_list",  access the P3 random forest ROC object as follows: 
# model_results_list$rf$P3$roc
###############################################################################
process_models <- function(train_test_list){
  # extract the individual train and test frames out of the master list
  class_names <- train_test_list[["class_names"]]
  train_set <- train_test_list[["train_set"]]
  train_class <- train_test_list[["train_class"]]
  test_set <- train_test_list[["test_set"]]
  test_class <- train_test_list[["test_class"]]
  train_set_rf <- train_test_list[["train_set_rf"]]
  train_class_rf <- train_test_list[["train_class_rf"]]
  indep_set <- train_test_list[["indep_set"]]
  
  # create empty lists for each model type
  knn_list <- list()
  rf_list <- list()
  
  # for each phenotype class in the class list
  for(phen_class in class_names){
    # extract the split training class vector from the training class frame
    train_class_n <- select_(train_class, .dots = phen_class)
    train_class_n <- train_class_n[[1]]
    
    # extract the random forest class vector from the RF class frame
    train_class_rf_n <- select_(train_class_rf, .dots = phen_class)
    train_class_rf_n <- train_class_rf_n[[1]]
    
    # if knn model processing is enabled
    if (knn_enabled){
      cat("Running KNN Model for class", phen_class, "\n")
      
      # run the KNN model
      start_time <- Sys.time()
      knn_results_list <- 
        process_knn(train_set, train_class_n, test_set)
      duration <- difftime(Sys.time(), start_time, units = "mins")
      cat("Run Time:", round(as.numeric(duration), 2), "minutes\n")
      
      # save the model results in the KNN model list
      knn_list[[phen_class]] <- knn_results_list
    }
    
    # if random forest processing is enabled
    if (rf_enabled){
      # save the start time
      start_time <- Sys.time()
      
      # run the random forest model
      cat("Running Random Forest Model for class", phen_class, 
          "-", format(start_time), "\n")
      rf_results_list <- process_rf(train_set_rf, train_class_rf_n, 
                                    indep_set)
      duration <- difftime(Sys.time(), start_time, units = "mins")
      cat("Run Time:", round(as.numeric(duration), 2), "minutes\n")
      
      # save the model results in the random forest model list
      rf_list[[phen_class]] <- rf_results_list  
    }
  }
  
  # concatenate the model lists together and return to the caller
  ret_model_results_list <- list()
  if (knn_enabled){
    ret_model_results_list$knn <- knn_list
  }
  
  if (rf_enabled){
    ret_model_results_list$rf <- rf_list
    ret_model_results_list$rf_drug_id <- train_test_list$train_set_drug_id_rf
  }
  
  return(ret_model_results_list)
}

###############################################################################
# Function Name: get_auc_results
# Description: Helper function which extracts an AUC results vector from
#              the model results list.
# Inputs: 
#   model_results_list: The master model results list.
# Returns a vector containing the AUC scores for each phenotype.
#
###############################################################################
get_auc_results <- function(model_results_list){
  auc_scores <- lapply(model_results_list$rf, function(i){
    i[["auc_score"]]
  })
  auc_scores <- unlist(auc_scores)
  auc_scores_df <- 
    data.frame(phenotype = names(auc_scores), auc_scores = auc_scores)
  rownames(auc_scores_df) <- NULL
  return(auc_scores_df)
}

###############################################################################
# Function Name: get_roc_results
# Description: Helper function which extracts an ROC results vector from
#              the model results list.
# Inputs: 
#   model_results_list: The master model results list.
# Returns a vector containing the ROC prediction objects for each phenotype.
#
###############################################################################
get_roc_results <- function(model_results_list){
  # extract a list of ROC objects from the model results list
  roc_results <- lapply(model_results_list$rf, function(i){
   i[["train_predict"]]
  })
  roc_results <- unlist(roc_results)
  return(roc_results)
}

###############################################################################
# Function Name: get_indep_predictions
# Description: Helper function which extracts the independent data set 
#              predictions from the model results.
# Inputs: 
#   model_results_list: The master model results list.
#   drug_ids: A vector containing the independent set drug IDs
#
# Returns a data frame containing the independent data set predictions.
#
###############################################################################
get_indep_predictions <- function(model_results_list, drug_ids){
  # get a list of prediction objects
  indep_predictions <- lapply(model_results_list$rf, function(i){
    i[["indep_predict"]]
  })
  
  # convert the prediction object list to a list of data frames.
  indep_predictions_list <- lapply(indep_predictions, function(i){
    data.frame(unlist(i))
  })
  
  # convert the list of data frames to a single data frame.
  indep_predictions_df <- data.frame(indep_predictions_list)
  
  # add the drug ID column to the data frame.
  indep_predictions_df <- indep_predictions_df %>%
    mutate(ID = drug_ids) %>%
    select(ID, everything())
  
  # return the data frame to the caller.
  return(indep_predictions_df)
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

# read the training data set
yactives_train <- read_csv("yactives_training_set_2.csv")

# read the independent test data
yactives_indep <- read_csv("CTRPv2_test_set.csv")

# get the training and test sets
train_prop <- 0.8
train_test_list <- 
  create_train_test_indep_sets(yactives_train, train_prop, yactives_indep)

# Proecess each model
# Notes on how to use the model_results_list.  The model results 
# list is a list of lists of lists. The top level is a list of models 
# like knn and rf.  The second level is a list of phenotypes like P1, P2, 
# P3, ... .  The third level is a list of model summary parameters like 
# fit object, predict object, and AUC score.  To access model parameters, 
# one can use the '$' notation.  For example, access the random forest 
# P3 ROC object as follows: model_results_list$rf$P3$roc
model_results_list <- process_models(train_test_list)

############### Save the results #####################
# save the master model_results_list to a file
save(model_results_list, file = "model_results_list.Rdata")

# extract a vector of AUC scores from the model results list and 
# save the results to a file.
auc_scores <- get_auc_results(model_results_list)
save(auc_scores, file = "auc_scores.Rdata")

# extract a vector of ROC predictions from the model results list
# and save the results to a file.
roc_results <- get_roc_results(model_results_list)
save(roc_results, file = "roc_results.Rdata")

# extract the independent test set predictions and save to a csv file
indep_predictions_df <- 
  get_indep_predictions(model_results_list, yactives_indep$ID)
write.csv(indep_predictions_df, file = "predictions.csv")
