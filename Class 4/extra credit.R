# lab 4 extra credit
source("http://thegrantlab.org/misc/cdc.R")
tail(cdc$height,20)
plot(cdc$height,cdc$weight)
weight_kg <- cdc$weight * 0.454
bmi <- (weight_kg)/(height_m^2)
table(cdc[bmi >= 30,9])
