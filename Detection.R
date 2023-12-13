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

#Load the the anomalous data sets
anomdata1 <- read.csv("Datasets/Anamalous_Datasets/DataWithAnomalies1.txt")
anomdata2 <- read.csv("Datasets/Anamalous_Datasets/DataWithAnomalies2.txt")
anomdata3 <- read.csv("Datasets/Anamalous_Datasets/DataWithAnomalies3.txt")

#Linearly interpolate all the missing values
anomdata1 <- na.locf(anomdata1)
anomdata1$Day_of_Week <-wday(dmy(anomdata1$Date))
anomdata1$Time_in_minutes <- as.ITime(anomdata1$Time)

target1 <- anomdata1[anomdata1$Day_of_Week == 4 &
                       anomdata1$Time_in_minutes >= as.ITime("13:00:00") &
                       anomdata1$Time_in_minutes < as.ITime("17:00:00"),
                     c("Date","Global_intensity","Global_active_power","Time_in_minutes")
]

time_stamps_count <- table(target1$Time_in_minutes)
model1 <- depmix(list(Global_intensity~1,Global_active_power~1), data=target1,
                nstates=18,
                family=list(gaussian(),gaussian()),
                ntimes=rep(240,52)
                )
fm1 <- fit(model1)
logLik(fm1)   #-25188.51 (18 states: -21138.24)


#Anomalous Data Set2
anomdata2 <- na.locf(anomdata2)
anomdata2$Day_of_Week <-wday(dmy(anomdata2$Date))
anomdata2$Time_in_minutes <- as.ITime(anomdata2$Time)

target2 <- anomdata2[anomdata2$Day_of_Week == 4 &
                       anomdata2$Time_in_minutes >= as.ITime("13:00:00") &
                       anomdata2$Time_in_minutes < as.ITime("17:00:00"),
                     c("Date","Global_intensity","Global_active_power","Time_in_minutes")
]

time_stamps_count <- table(target2$Time_in_minutes)
model2 <- depmix(list(Global_intensity~1,Global_active_power~1), data=target2,
                 nstates=18,
                 family=list(gaussian(),gaussian()),
                 ntimes=rep(240,52)
)
fm2 <- fit(model2)
logLik(fm2)     #-25626.92 (with 18 states: -21823.92)



#Anomalous Data Set3
anomdata3 <- na.locf(anomdata3)
anomdata3$Day_of_Week <-wday(dmy(anomdata3$Date))
anomdata3$Time_in_minutes <- as.ITime(anomdata3$Time)

target3 <- anomdata3[anomdata3$Day_of_Week == 4 &
                       anomdata3$Time_in_minutes >= as.ITime("13:00:00") &
                       anomdata3$Time_in_minutes < as.ITime("17:00:00"),
                     c("Date","Global_intensity","Global_active_power","Time_in_minutes")
]

time_stamps_count <- table(target3$Time_in_minutes)
model3 <- depmix(list(Global_intensity~1,Global_active_power~1), data=target3,
                 nstates=18,
                 family=list(gaussian(),gaussian()),
                 ntimes=rep(240,52)
)
fm3 <- fit(model3)
logLik(fm3)     #-54952.81 (with 18 states: -52231.92)
