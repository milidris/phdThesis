###############################################
# Packages import
library(keras)

###############
# Data Import

load(file="useCases/acousticFireExtinguisher/data/AFE_data.RData")

##################
# Data Preparation

# One-Hot the categorical variable
data.AFE<-fire[,-2]
data.AFE$Gasoline<-as.numeric(fire$Fuel==1)
data.AFE$Kerosene<-as.numeric(fire$Fuel==2)

# Split Train/Test
set.seed(12345)
n.val = round(nrow(data.AFE)*0.05) #95% of data for training
idx.val=sample(seq(1, nrow(data.AFE)), n.val)
idx.train=setdiff(seq(1, nrow(data.AFE)), idx.val)

# Setup of keras-compliant data structures
X.train<-as.matrix(data.AFE[idx.train,-6])
Y.train<-to_categorical(as.matrix(data.AFE[idx.train,6]))

X.val<-as.matrix(data.AFE[idx.val, -6])
Y.val<-to_categorical(as.matrix(data.AFE[idx.val, 6]))

#########################
# NN-model architecture

input <- layer_input(shape=7) #7 features as inputs
output <- input %>%
  layer_dense(units = 100, activation = "relu") %>% # One layer
  layer_dense(units = 2, activation = 'softmax')

model.AFE <- keras_model(input, output)

model.AFE %>%
  compile(
    loss = "binary_crossentropy",
    optimizer = optimizer_adam(lr = 0.0001), #Learning rate setup
    metrics = c("accuracy")
  )

early_stop <- callback_early_stopping(monitor = "val_loss", patience = 50)
# Stop the training if no improvement on the validation loss for 50 epochs

# Model training
set.seed(54321)
model.AFE %>% fit(
  x = X.train,
  y = Y.train,
  batch_size=64,
  epochs = 500,
  verbose = 2,
  validation_data=list(X.val, Y.val),
  callbacks=list(early_stop)
)

## 94.53% Train accuracy
## 93.51% Test accuracy

# Saving the model
save_model_hdf5(model.AFE, "useCases/acousticFireExtinguisher/model/model_AFE.hdf5")

save(data.AFE, X.train, X.val, Y.train, Y.val, file="useCases/acousticFireExtinguisher/model/AFE_trainingData.RData")
