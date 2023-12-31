###########################
# Libraries
library(readxl)

# Extract data from initial .xlsx file
fire <- read_excel("useCases/acousticFireExtinguisher/data/fireData.xlsx", sheet = "fire")

###########################
# Data Preparation
fire<-fire[fire$FUEL!="lpg",] #We don't keep "lpg" fuel
fire$SIZE=as.numeric(fire$SIZE)
fire$SIZE[fire$SIZE==1]<-7
fire$SIZE[fire$SIZE==2]<-12
fire$SIZE[fire$SIZE==3]<-14
fire$SIZE[fire$SIZE==4]<-16
fire$SIZE[fire$SIZE==5]<-20

colnames(fire)=c("TankSize", "Fuel","Distance" ,"Decibel", "Airflow", "Frequency", "Y")

fire$Fuel<-factor(fire$Fuel,
                  levels=list("gasoline",
                              "kerosene",
                              "thinner"),
                  labels=c(1,2,3)
)


###########################
# Saving Data
save(fire, file="useCases/acousticFireExtinguisher/data/AFE_data.RData")
