#For reproducible results for the model
set.seed(1)

# List of required packages
packages <- c("ggplot2","lubridate","data.table","depmixS4","dplyr","zoo")

# Check, install required packages
for (package in packages) {
  if (!require(package, character.only = TRUE)) { #character.only is used to ensure we check for the package name as a string
    install.packages(package)
  }
  library(package, character.only = TRUE)
}
# Load the data and extract the non-NA values from it
termProjData <- read.csv("Datasets/TermProjectData.txt", header = TRUE)

#na.locf does not fill in the first missing value
termProjData["Global_active_power"][1,] <- 5.36000

#Apply linear interpolation to all columns with missing values
interpolatedData <- na.locf(termProjData)

#Add a column to indicate the day of the week
interpolatedData$Day_of_week <- wday(dmy(interpolatedData$Date))

#Add a column for time in minutes(integer object)
interpolatedData$Time_in_minutes <- as.ITime(interpolatedData$Time)

#Change the date to date object in interpolated data
interpolatedData$Date <- as.Date(interpolatedData$Date, format = "%d/%m/%Y")

# Get the PCA of the 7 responses with scaled data
pca <- prcomp(interpolatedData[, 3:9], scale = TRUE)

# Put the data in a nice format for ggplot to use

pca.data <- data.frame(Sample = seq_len(nrow(pca$x)),
                       X = pca$x[, 1],
                       Y = pca$x[, 2])
# Calculate the percentage of each PC
pca.var <- pca$sdev^2

pca.var.per <- round(pca.var/sum(pca.var)*100, 1)
barplot(pca.var.per, main = "Scree plot", xlab = "Principal Component",
        ylab = "Percent Variation")

# Plot the PCA graph for the 7 responses
ggplot(data = pca.data, aes(x = X, y = Y, label = Sample)) +
  geom_point() +
  xlab(paste0("PC1 - ", pca.var.per[1], "%")) +
  ylab(paste0("PC2 - ", pca.var.per[2], "%")) +
  theme_bw() +
  ggtitle("PCA Graph")

# Plot loading scores on PCA graph
biplot(pca)


# Get the info for which vars push the data to the left and which vars push the data to the right
#Information on which variables create the two PCA clusters
loadingScores <- pca$rotation[, 1]
varScores <- abs(loadingScores)
varScoresRanked <- sort(varScores, decreasing = TRUE)
top7Vars <- names(varScoresRanked[1:7])
pca$rotation[top7Vars, 1] 

#Based on the loading scores of the PCA we select Global_intensity and Global_active_power
#Select a training data i.e. the first two years
#Day selected: Wednesday.  4 hours time window from 1:00pm to 5:00pm

target_data <- interpolatedData[interpolatedData$Day_of_week==4 &
                                  interpolatedData$Time_in_minutes >= as.ITime("13:00:00") &
                                  interpolatedData$Time_in_minutes < as.ITime("17:00:00"),
                                c("Date","Global_intensity","Global_active_power","Time_in_minutes"),
]

training_data <- target_data[target_data$Date < as.Date("2009-01-01"), ]

training_data[2:3] <- scale(training_data[2:3])[,1:2]

#Identify the distribution of the chosen variables
#Plot Global_intensity against time
ggplot(data = training_data, aes(x = Time_in_minutes, y = Global_intensity))+
  geom_point(size = 0.5) +#Use points to represent data
  labs(title = "Distribution of Global_intensity over Time",
       x = "Time in minutes",
       y = "Glabal_intensity") +
  theme_minimal()

#Plot the Global_active_power against time
ggplot(data = training_data, aes(x = Time_in_minutes, y = Global_active_power))+
  geom_point(size = 0.5) +#Use points to represent data
  labs(title = "Distribution of Global_active_power over Time",
       x = "Time in minutes",
       y = "Glabal_active_power") +
  theme_minimal()

#Group the data points that occur at the same time
training_data <- training_data[order(training_data$Time_in_minutes), ]

# test run to see if n times works
# model <- depmix(list(Global_intensity~1,Global_active_power~1), data=training_data,
#                 nstates=4,
#                 family=list(gaussian(),gaussian()),
#                 ntimes=rep(240,107)
# )

#Define a function to fit the model given different states
fit_model <- function(n_states){
  model <- depmix(list(Global_intensity~1,Global_active_power~1), 
                  data=training_data,nstates=n_states, 
                  family=list(gaussian(),gaussian()),
                  ntimes=rep(240,107)
  )
  fit <- fit(model)
  
  #Return a dataframe with the BIC, number of states and logliklihood
  data.frame(
    n_states = n_states,
    BIC = BIC(fit),
    log_likelihood = logLik(fit)
  ) 
}

#Set the range of states to experiment in i.e 4 to 18
state_range <- 4:18 

results <- data.frame(n_states= integer(),
                      BIC = numeric(),
                      log_likeihood = numeric())

#Loop through different states and store the results in the dataframe
for (n_states in state_range) {
  results <- rbind(results, fit_model(n_states))
}

#Plot BIC versus likelihood to find the optimal number of states for the HMM
ggplot(data = results, aes(x=BIC, y = log_likelihood))+
  geom_text(label = results$n_states, size = 3) +
  labs(title = "Loglikelihood versus BIC",
       x = "BIC",
       y = "Logliklihood")+
  theme_minimal()

#From the plot we see that the most optimal states are 9,10 & 11

#Selecting test data
test_data <- target_data[target_data$Date >= as.Date("2009-01-01"), ]

#Find the occurrence of each time stamp in test_data
time_stamps_count <- table(test_data$Time_in_minutes)
#print(time_stamps_count)

loglikelihoods <- numeric()
#Compute the likelihood for using HMM with states n = 9,10 & 11
for (i in 16) {
  model <- depmix(list(Global_intensity~1,Global_active_power~1), data=test_data,
                  nstates=i,
                  family=list(gaussian(),gaussian()),
                  ntimes=rep(240,47)
  )
  fm <- fit(model)
  
  loglikelihoods <- c(loglikelihoods,logLik(fm))
}

#From the test data we get
# { 9:-19797.73 , 10:-20209.08, 11:-18780.91}
# { 18: -15201.49 }
#No overfitting, higher loglikelihood than the trainign data

#The ideal number of states from the training data n = 11


#Using the model on the anomalous data sets


ggplot(data = results, aes(x = BIC, y = log_likelihood))+
  geom_text(label = results$n_states, size = 3) +
  labs(title = "Loglikelihood versus BIC",
       x = "BIC",
       y = "Logliklihood")+
  theme_minimal()

ggplot(data = results, aes(x = n_states, y = log_likelihood))+
  geom_point() +
  labs(title = "Loglikelihood versus N States",
       x = "BIC",
       y = "Logliklihood")+
  theme_minimal()

ggplot(data = results, aes(x = n_states, y = BIC))+
  geom_point() +
  labs(title = "BIC versus N States",
       x = "BIC",
       y = "Logliklihood")+
  theme_minimal()

# Get normalized loglikelihood data
results$log_likelihood/(107*240)

# Normalized training loglikelihood
# 15: -1.0213514
# 16: -1.0131124
# 18: -0.9992406

# Normalized test loglikelihood
# 15: -1.481501
# 16: -1.522597
# 18: -1.34765