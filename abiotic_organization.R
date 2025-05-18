########################################################
#DC | 240224
########################################################
library(readxl)
library(plyr)
library(readr)
library(data.table)
library(dplyr)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(ggpubr)
library(gridExtra)
library(devtools)
library(plotrix)
library(lme4)
library(gplots) # This package provides the heatmap function
library(ggrepel)

#######################
#Bringing the abiotic and ecosystem functioning data in one script comparable to the trait data
#######################
setwd("/Users/damlacinoglu/sarah trait paper check")

merged = read.table(file = "MERGED_231122")
head(merged)

abiotic_merged = data.frame(matrix(ncol = 15, nrow = nrow(merged)))
abiotic_merged[1:nrow(abiotic_merged),1:5] = merged[1:nrow(merged),1:5]
colnames(abiotic_merged)[1:5] = colnames(merged)[1:5]
colnames(abiotic_merged)[6:15]= c("temp", "moist","total_cover","total_biomass","sp_cover","sp_biomass","light_ground1","light_ground2","light_canopy","nh4")

temperature = data.frame(read_excel("TEMPERATURE_END_2022.xlsx"))
moisture = data.frame(read_excel("MOISTURE_END_2022.xlsx"))
percent_cover = data.frame(read_excel("COVER_END_2022.xlsx"))
biomass = data.frame(read_excel("BIOMASS_2022.xlsx"))
light = data.frame(read_excel("LIGHT_2023.xlsx"))
soil_nitrogen = data.frame(read.csv("soilNH4240319_SKO.csv"))

length(unique(abiotic_merged$plot))
unique(temperature$plot) %in% unique(percent_cover$Plot)
unique(percent_cover$Plot) %in% unique(light$Plot.ID)
unique(temperature$plot) %in% unique(light$Plot.ID)
unique(abiotic_merged$plot) %in% unique(light$Plot.ID)
unique(abiotic_merged$plot) %in% unique(percent_cover$Plot)
unique(abiotic_merged$plot) %in% unique(temperature$plot)

abiotic_merged$plot = gsub("sp","",abiotic_merged$plot )

for(i in 1:nrow(abiotic_merged)) {
	plot = abiotic_merged[i,"plot"]
	species = abiotic_merged[i,"species"]
	
	 abiotic_merged[i,"temp"] = as.numeric(sum(temperature[temperature[,"plot"] == plot,3:5])/3)
	 abiotic_merged[i,"moist"] = as.numeric(moisture[moisture[,"plot"] == plot,4])
	 
	 	 #total cover is total - weeds (total and weeds exist for all plots)
	abiotic_merged[i,"total_cover"] = percent_cover[percent_cover[,"Plot"] == plot & percent_cover[,"Species.Code"] == ".total","Out.net..."] 
	 abiotic_merged[i,"sp_cover"] = percent_cover[percent_cover[,"Plot"] == plot & percent_cover[,"Species.Code"] == species,"Out.net..."]
	 	  
	 #total biomass is total + weeds for monocultures and total (weed already included in 12sp)
	    #for anything that has a bag 2.
	    abiotic_merged[i,"total_biomass"] = sum(na.omit(as.numeric(biomass[biomass[,"Plot.name"] == plot & biomass[,"Quadrat"] == "Open","Total.weight..g."])))
		 abiotic_merged[i,"sp_cover"] = percent_cover[percent_cover[,"Plot"] == plot & percent_cover[,"Species.Code"] == species,"Out.net..."] / 100
	 abiotic_merged[i,"sp_biomass"] = abiotic_merged[i,"total_biomass"] * abiotic_merged[i,"sp_cover"]
	  
 abiotic_merged[i,"light_ground1"] = as.numeric(light[light[,"Plot.ID"] == plot & light[,"Quadrat"] == "OPEN","Ground.1"])
  abiotic_merged[i,"light_ground2"] = as.numeric(light[light[,"Plot.ID"] == plot & light[,"Quadrat"] == "OPEN","Ground.2"])    
  abiotic_merged[i,"light_canopy"] = as.numeric(light[light[,"Plot.ID"] == plot & light[,"Quadrat"] == "OPEN","Above.Canopy"])
     abiotic_merged[i,"nh4"] = as.numeric(soil_nitrogen[soil_nitrogen[,"plot"] == plot,"nh4.mean"])   
}

write.table(abiotic_merged, file = "ABIOTIC_250304")

nrow(abiotic_merged)
head(abiotic_merged)
tail(abiotic_merged)



# Summarizing columns from 6th onward by 'plot'
abiotic_summary <- abiotic_merged %>%
  group_by(plot) %>%
  summarize(across(temp:nh4, mean, na.rm = TRUE))

# View the summarized data
print(abiotic_summary, n = Inf)






