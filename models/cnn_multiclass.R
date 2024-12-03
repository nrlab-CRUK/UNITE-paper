
library(tidyverse)
library(keras)
library(tensorflow)
library(CatEncoders)
library(superml)



plot_path <- "/mnt/scratchc/nrlab/wang04/ulyses"
plot_file_name <-file.path(plot_path, "cnn_training_history_multi_class.pdf")
cohorts <- c("Breast", "Prostate", "Lung", "Healthy")

tensors <- readRDS("/mnt/scratchc/nrlab/wang04/ulyses/tensors.rds")
labels <- readRDS("/mnt/scratchc/nrlab/wang04/ulyses/labels.rds")




# functions ---------------------------------------------------------------

class_weights <- function(labels, case = 1, mu = 0.15) {
  # do bin count
  weights <- table(labels)
  if (case == 1){
    # sklearn.utils.class_weight.compute_class_weight approach
    # http://scikit-learn.org/stable/modules/generated/sklearn.utils.class_weight.compute_class_weight.html
    weights <- sum(weights) / (length(weights) * weights)
  } else if (case == 2) {
    weights <- log(mu * sum(weights) / weights)
    weights <- ifelse(weights < 1, 1, weights)
  } else if (case == 3) {
    weights <- ceiling(max(weights) / weights)
  } else {
    weights <- weights / sum(weights)
  }
  # create and return list
  setNames(as.list(weights), names(weights))
}

# pre-processing x and y --------------------------------------------------

# subset input based on cohorts
ids <- which(unlist(labels) %in% cohorts)
tensors <- tensors[ids]
labels <- labels[ids]

# shuffle input
tmp <- sample(1:length(labels))
tensors <- tensors[tmp]
labels <- labels[tmp]

# make x and y
train_x <- abind::abind(tensors, along = 0)
lbl <- superml::LabelEncoder$new()
train_y <- lbl$fit_transform( unlist(labels)) %>%
  keras::to_categorical(num_classes = lbl$encodings %>% length())

# set class weight
class_weight <- class_weights(labels = lbl$fit_transform(unlist(labels)))
input_shape <- dim(tensors[[1]])
n_class <- lbl$encodings %>% length()

# convolutional layer
model <- keras_model_sequential() %>% 
  layer_conv_2d(filters = 32, kernel_size = c(3, 3), activation = "relu", 
                input_shape = input_shape) %>% 
  layer_max_pooling_2d(pool_size = c(3,2)) %>% 
  layer_conv_2d(filters = 32, kernel_size = c(3, 3), activation = "relu") %>% 
  layer_max_pooling_2d(pool_size = c(2,2)) %>%
  layer_conv_2d(filters = 32, kernel_size = c(3, 3), activation = "relu") %>% 
  layer_max_pooling_2d(pool_size = c(2,2)) %>%
  layer_conv_2d(filters = 16, kernel_size = c(3, 3), activation = "relu") %>% 
  layer_max_pooling_2d(pool_size = c(2,2)) 


model %>% 
  layer_flatten() %>%
  # layer_global_average_pooling_2d() %>% 
  layer_dense(units = 32, activation = "relu") %>%
  layer_dropout(rate = 0.2) %>%
  layer_dense(units = n_class, activation = "softmax")

model %>% compile(
  optimizer = optimizer_adam(learning_rate = 0.001),
  loss = "categorical_crossentropy",
  metrics = list("accuracy", 
                 "AUC", 
                 metric_sensitivity_at_specificity(specificity = 0.95), 
                 metric_precision(), 
                 metric_recall())
)

# train the model

history_cnn <- model %>% 
  fit(
    x = train_x, y = train_y,
    epochs = 40,
    batch_size = 50,
    verbose = 1,
    callbacks = callback_tensorboard(),
    validation_split = 0.2,
    shuffle = TRUE,
    class_weight = class_weight
  )


# save the model ----------------------------------------------------------

#plot the history
pdf(file = plot_file_name )
plot(history_cnn)
dev.off()

# save the r image --------------------------------------------------------
save.image(file = "cnn_multiclass.rda")


# transfer ----------------------------------------------------------------

model <- application_vgg19(weights = NULL, 
                           classes = n_class, 
                           include_top = TRUE, 
                           input_shape = input_shape )


model %>% compile(
  optimizer = "adam",
  loss = "categorical_crossentropy",
  metrics = list("accuracy", 
                 "AUC", 
                 metric_sensitivity_at_specificity(specificity = 0.95), 
                 metric_precision(), 
                 metric_recall())
)

history <- model %>% 
  fit(
    x = train_x, y = train_y,
    epochs = 40,
    batch_size = 50,
    verbose = 1,
    callbacks = callback_tensorboard(),
    shuffle = TRUE,
    class_weight = class_weight,
    validation_split = 0.2
    
  )


