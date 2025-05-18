########################################################
#DC | 231031
#Here, I'm trying to wrangle the SLA, Leaf N and C, and height data from 2022 into one dataset. This is all the data coming from the individuals. 
#In addition, I also separate out what's left for the height data in another dataset.
########################################################
install.packages("readxl")
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

setwd("/Users/damlacinoglu/Desktop/final_bfl_data")
sla = read.table("SLA_2022.txt", sep = ",")
nitrogen = data.frame(read_excel("LEAF_C_N_2022_2023.xlsx"))
height = data.frame(read_excel("HEIGHT_2022.xlsx"))
all_codes = data.frame(read_excel("ALL_CODES_2022.xlsx"))

head(sla); head(nitrogen); head(height); head(all_codes)

all_codes = cbind(all_codes, LA = rep(NA, nrow(all_codes)), LMA = rep(NA, nrow(all_codes)), d15N = rep(NA, nrow(all_codes)), d13C = rep(NA, nrow(all_codes)), percN = rep(NA, nrow(all_codes)), percC = rep(NA, nrow(all_codes)), veg_height = rep(NA, nrow(all_codes)), rep_height = rep(NA, nrow(all_codes)))

nitrogen$Sample.ID <- gsub("APRU", "ARPU", nitrogen$Sample.ID)
height$indvidual <- gsub("APRU", "ARPU", height$indvidual)

nitrogen_22 = nitrogen[1:418,]
nitrogen_23 =  nitrogen[419:nrow(nitrogen),]

########################################################
#Individual allocation data
########################################################

for(i in 1:nrow(all_codes)) {
	 ind = all_codes[i, "indvidual"]
	 
	 if(nrow(sla[sla[,"indvidual"] == ind,]) == 1) {
	 	  all_codes[i,"LA"] = sla[sla[,"indvidual"] == ind,"area_2"]
	      all_codes[i,"LMA"] = 1 / sla[sla[,"indvidual"] == ind,"sla"]
	 } else {
	 	  all_codes[i,"LA"] = NA; all_codes[i,"LMA"] = NA }

 	  if(nrow(nitrogen_22[nitrogen_22[,"Sample.ID"] == ind,]) == 1) {
 	  	 all_codes[i,"d15N"] = nitrogen_22[nitrogen_22[,"Sample.ID"] == ind,"d15N..permil."]
 	  	 all_codes[i,"d13C"] = nitrogen_22[nitrogen_22[,"Sample.ID"] == ind,"d13C..permil."]
	     all_codes[i,"percN"] = nitrogen_22[nitrogen_22[,"Sample.ID"] == ind,"N."]
	     all_codes[i,"percC"] = nitrogen_22[nitrogen_22[,"Sample.ID"] == ind,"C."]
 	 } 
 	 
 	  if(nrow(nitrogen_22[nitrogen_22[,"Sample.ID"] == ind,]) > 1) {
 	  	#some nitrogen might have multiple values, average them if they come from the same sample (duplicates and resends)
 	  	print(i)
 	  	print("*****")
 	  	 all_codes[i,"d15N"] = mean(nitrogen_22[nitrogen_22[,"Sample.ID"] == ind,"d15N..permil."])
 	  	 all_codes[i,"d13C"] = mean(nitrogen_22[nitrogen_22[,"Sample.ID"] == ind,"d13C..permil."])
	     all_codes[i,"percN"] = mean(nitrogen_22[nitrogen_22[,"Sample.ID"] == ind,"N."])
	     all_codes[i,"percC"] = mean(nitrogen_22[nitrogen_22[,"Sample.ID"] == ind,"C."])
 	 } 
 	 
 	 if(nrow(nitrogen_22[nitrogen_22[,"Sample.ID"] == ind,]) < 1) {
 	  	  all_codes[i,"d15N"] = NA; all_codes[i,"d13C"] = NA
	     all_codes[i,"percN"] = NA; all_codes[i,"percC"] = NA
	  } 
 	  
 	 if(nrow(height[height[,"indvidual"] == ind,]) == 1) {
 	  	 all_codes[i,"veg_height"] = height[height[,"indvidual"] == ind,"vegetative.height."]
 	  	 all_codes[i,"rep_height"] = height[height[,"indvidual"] == ind,"reproductive.height."]
 	 } else {
 	  	 all_codes[i,"veg_height"] = NA; all_codes[i,"rep_height"] = NA
 	  }
	
}

all_codes = cbind(all_codes, rep(2022, nrow(all_codes)))
colnames(all_codes)[14] = c("year")

#height should have left overs 
`%notin%` <- Negate(`%in%`)

leftover_height = height[height[,"indvidual"] %notin% all_codes[,"indvidual"],]
for(i in 1:nrow(leftover_height)) {
	all_codes = rbind(all_codes,c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA))
	all_codes[432+i,"indvidual"] = leftover_height[i,"indvidual"]
	all_codes[432+i,"veg_height"] = leftover_height[i,"vegetative.height."]
		 all_codes[432+i,"rep_height"] = leftover_height[i,"reproductive.height."]
 all_codes[432+i,"year"] = 2022	     

}

nrow(all_codes)

#Add a marker to mark 2023 nitrogen stuff to remove from individual allocation
colnames(all_codes)[14] = c("year")
for(i in 1:nrow(nitrogen_23)) {
	all_codes = rbind(all_codes,c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA))
	all_codes[1806+i,"indvidual"] = nitrogen_23[i,2]
		 all_codes[1806 +i,"d15N"] = nitrogen_23[i,"d15N..permil."]
 	  	 all_codes[1806 +i,"d13C"] = nitrogen_23[i,"d13C..permil."]
	     all_codes[1806 +i,"percN"] = nitrogen_23[i,"N."]
	     all_codes[1806 +i,"percC"] = nitrogen_23[i,"C."]
 all_codes[1806 +i,"year"] = 2023	     
}

#Add other stuff's 2023 data - can't do for height but for SLA? ADD 2023 SLA LATER!!!
list_to_add = all_codes[all_codes[,"year"] == 2022 & is.na(all_codes[,"LMA"]),]

#Add the info from the first couple of columns back in for the nitrogen and height stuff
for(i in 433:nrow(all_codes)) {
	code = all_codes[i,1]
	if(grepl("W", code)) {all_codes[i,"treatment"] = "watered"} else {all_codes[i,"treatment"] = "control"}
	if(grepl("mon", code)) {all_codes[i,"richness"] = 1} else {all_codes[i,"richness"] = 12}
	all_codes[i,2] = strsplit(code, "[-]")[[1]][1]
	if(all_codes [i,"treatment"] == "control") {	all_codes[i,3] =  strsplit(code, "[-]")[[1]][2]} else {all_codes[i,3] = paste(strsplit(code, "[-]")[[1]][2],"-W", sep = '')}
}

for(i in 433:nrow(all_codes)) {
		code = all_codes[i,3] 

		if(all_codes[i,"richness"] == 12) {
					#Add sp for 12 sp plot codes
			if(grepl("12", code)) {all_codes[i,3] = gsub("12","12sp", code)}} else { 
							#Capitalize monoculture names
		#Add unique codes for monoculture plots
		if(all_codes[i,2] == "ARPU") {all_codes[i,3] = gsub("mon","MON1", code)}
		if(all_codes[i,2] == "LUTE") {all_codes[i,3] = gsub("mon","MON5", code)}
		if(all_codes[i,2] == "GAPU") {all_codes[i,3] = gsub("mon","MON9", code)}
		if(all_codes[i,2] == "MOCI") {all_codes[i,3] = gsub("mon","MON8", code)}
		if(all_codes[i,2] == "MOPU") {all_codes[i,3] = gsub("mon","MON7", code)}
		if(all_codes[i,2] == "COTI") {all_codes[i,3] = gsub("mon","MON10", code)}
		if(all_codes[i,2] == "PHCO") {all_codes[i,3] = gsub("mon","MON4", code)}
		if(all_codes[i,2] == "CEAM") {all_codes[i,3] = gsub("mon","MON11", code)}
		if(all_codes[i,2] == "BOCU") {all_codes[i,3] = gsub("mon","MON2", code)}
		if(all_codes[i,2] == "BOGR") {all_codes[i,3] = gsub("mon","MON3", code)}
		if(all_codes[i,2] == "DEIL") {all_codes[i,3] = gsub("mon","MON6", code)}
		if(all_codes[i,2] == "IPRU") {all_codes[i,3] = gsub("mon","MON12", code)}}
}

unique(all_codes[,3])
all_codes[all_codes[,3] == "NA",]
all_codes = all_codes[!is.na(all_codes[,3]),]

write.table(all_codes, file = "ABOVEGROUND_2022")

########################################################
#Merging above and below ground data
########################################################

setwd("/Users/damlacinoglu/Desktop/final_bfl_data")
above = read.table("ABOVEGROUND_2022", sep = " ")
below = read.csv("raw_root_data_2022_sko.csv", sep = ",")

#have to edit the aboveground data to make it match below
	#change all dashes to underscores
	#make all Mons MONS
colnames(above)[1] = c("individual")
for(i in 1:nrow(above)) {
	name = above[i,"individual"]
	above[i,"individual"] = gsub(pattern = "mon", replacement = "MON", name)
	name = above[i,"individual"]
	above[i,"individual"] = gsub(pattern = "-", replacement = "_", name)
	if(above[i,"treatment"] == "control") {above[i,"treatment"] = "C"}
	if(above[i,"treatment"] == "watered") {above[i,"treatment"] = "W"}
}

write.table(above, file = "ABOVEGROUND_2022_v2")

merged = cbind(above, rep(NA, nrow(above)), rep(NA, nrow(above)), rep(NA, nrow(above)), rep(NA, nrow(above)), rep(NA, nrow(above)), rep(NA, nrow(above)), rep(NA, nrow(above)), rep(NA, nrow(above)), rep(NA, nrow(above)), rep(NA, nrow(above)), rep(NA, nrow(above)), rep(NA, nrow(above)), rep(NA, nrow(above)), rep(NA, nrow(above)), rep(NA, nrow(above)), rep(NA, nrow(above)), rep(NA, nrow(above)), rep(NA, nrow(above)), rep(NA, nrow(above)))
colnames(merged)[15:33] = c("family","growth.form","life.hist","Total.Root.Length.mm","Depth.mm","depth.cm","Average.Diameter.mm","fine.root.length.mm","fine.root.length.m","weight.srl","srl","shoot.weight","coarse.root.weight","fine.root.weight","total.root.weight","total.biomass","root.shoot.ratio","root.tissue.density","root.mass.fraction")

for(i in 1:nrow(merged)) {
	sp = merged[i,"species"]
		sp_data = below[below[,"species"] == sp,]
		merged[i,"family"] = sp_data[1,"family"]
		merged[i,"growth.form"] = sp_data[1,"growth.form"]
		merged[i,"life.hist"] = sp_data[1,"life.hist"]
	ind = merged[i,"individual"]
	if(ind %in% below[,"individual"]) {
		ind_data = below[below[,"individual"] == ind, ]
		merged[i,18:33] = ind_data[,9:24]
	}
}

write.table(merged, file = "MERGED_231122")

########################################################
#Plotting
########################################################
#check units and values of SLA. CHECK
#check height averages with the wildflower data

##### Traits by treatment graphs ####

#######################
#LMA
#######################

lma.wat <- merged_2022 %>%
drop_na(LMA) %>% 
  group_by(treatment) %>%
  summarise_at(vars(LMA), list(sd = sd,
                               se = std.error,  
                               LMA = mean))
lma.wat

lma.wat.sp <- merged_2022 %>%
drop_na(LMA) %>% 
  group_by(species, treatment) %>%
  summarise_at(vars(LMA), list(sd = sd,
                               se = std.error,  
                               LMA = mean))
lma.wat.sp


lma.g1.wat <- ggplot(data = lma.wat.sp,
                     aes(y = LMA, 
                         x = treatment,
                         col = treatment)) + 
  geom_point(size = 3,
             alpha = 0.2) +
  geom_line(aes(group = species),
            col = 'gray88') + 
  geom_point(data = lma.wat,
             aes(y = LMA,
                 x = treatment,
                 col = treatment),
             size = 2) +
  geom_errorbar(data = lma.wat,
                aes(ymin = LMA-se,
                    ymax = LMA+se),
                size = 1,
                width = 0.2) + 
  ylab('Leaf Mass per Area (kg/m^2)') +  
  xlab('Water Treatment') + 
  theme(text = element_text(size = 5)) + 
  theme_bw(base_size = 12) + 
  theme(axis.text.x = element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = 'bottom',
        strip.text = element_text(size = 10),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        legend.background = element_rect(fill='transparent'),
        legend.box.background = element_rect(fill='transparent')) + 
  scale_color_manual(values= c('skyblue2', 'skyblue4'),
                     labels=c('Ambient', 'Watered')) +
  labs(color = " ")
lma.g1.wat

#anova
lma.wat.aov <- aov(LMA~species*treatment, data = merged_2022)
summary(lma.wat.aov) # residuals > sub sq --> intra > inter
# species signicant, but not treatment and sp*treat


#######################
#VEG HEIGHT
#######################

merged_2022$veg_height = as.numeric(merged_2022$veg_height)
merged_2022$rep_height = as.numeric(merged_2022$rep_height)

veg.wat <- merged_2022 %>%
  drop_na(veg_height) %>% 
  group_by(treatment) %>%
  summarise_at(vars(veg_height), list(sd = sd,
                               se = std.error,  
                               height = mean))
veg.wat

veg.wat.sp <- merged_2022 %>%
  drop_na(veg_height) %>% 
  group_by(species, treatment) %>%
  summarise_at(vars(veg_height), list(sd = sd,
                               se = std.error,  
                               height = mean))
veg.wat.sp

rep.wat.sp <- merged_2022 %>%
  drop_na(rep_height) %>% 
  group_by(species, treatment) %>%
  summarise_at(vars(rep_height), list(sd = sd,
                               se = std.error,  
                               height = mean))
rep.wat.sp

veg.g1.wat <- ggplot(data = veg.wat.sp,
                     aes(y = height, 
                         x = treatment,
                         col = treatment)) + 
  geom_point(size = 3,
             alpha = 0.2) +
  geom_line(aes(group = species),
            col = 'gray88') + 
  geom_point(data = veg.wat,
             aes(y = height,
                 x = treatment,
                 col = treatment),
             size = 2) +
  geom_errorbar(data = veg.wat,
                aes(ymin = height-se,
                    ymax = height +se),
                size = 1,
                width = 0.2) + 
  ylab('Vegetative Height (cm)') +  
  xlab('Water Treatment') + 
  theme(text = element_text(size = 5)) + 
  theme_bw(base_size = 12) + 
  theme(axis.text.x = element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = 'bottom',
        strip.text = element_text(size = 10),
        panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        legend.background = element_rect(fill='transparent'),
        legend.box.background = element_rect(fill='transparent')) + 
  scale_color_manual(values= c('skyblue2', 'skyblue4'),
                     labels=c('Ambient', 'Watered')) +
  labs(color = " ")
veg.g1.wat

#anova
veg.wat.aov <- aov(height~species*treatment, data = merged_2022)
summary(veg.wat.aov) # residuals > sub sq --> intra > inter
# species signicant, but not treatment and sp*treat

#######################
#PERC CARBON
#######################


#######################
#PERC NITROGEN
#######################

#######################
#CORRELATIONS + HEAT MAPS
#######################
merged = read.table(file = "MERGED_231122")
head(merged)

install.packages("gplots")
library(gplots) # This package provides the heatmap function

# Set heatmap colors (optional)
my_palette <- colorRampPalette(c("blue", "white", "red"))(100)

#AMBIENT
cor_data = merged[1:432,]
cor_data = cor_data[cor_data[,"treatment"] == "C",]
cor_data_2 = cbind(cor_data[6:13], cor_data[18:33])

numeric_cols <- sapply(cor_data_2, is.numeric)
cor_data_2[, !numeric_cols] <- lapply(cor_data_2[, !numeric_cols], as.numeric)
#cor_data_2[,1:ncol(cor_data_2)] = scale(cor_data_2[,1:ncol(cor_data_2)])
cor_data_2[,1:ncol(cor_data_2)] = log(cor_data_2[,1:ncol(cor_data_2)])

correlation_matrix <- cor(cor_data_2, method = "pearson",use = "pairwise.complete.obs")

# Plotting the heatmap
quartz()
heatmap(
  correlation_matrix,
  col = my_palette,
  symm = TRUE, # Display the heatmap symmetrically
  margins = c(10, 10), # Add margins around the heatmap
  main = "Ambient",
   Rowv = NA, 
  Colv = NA
)

#WATER
cor_data = merged[1:432,]
cor_data = cor_data[cor_data[,"treatment"] == "W",]
cor_data_2 = cbind(cor_data[6:13], cor_data[18:33])

numeric_cols <- sapply(cor_data_2, is.numeric)
cor_data_2[, !numeric_cols] <- lapply(cor_data_2[, !numeric_cols], as.numeric)
cor_data_2[,1:ncol(cor_data_2)] = scale(cor_data_2[,1:ncol(cor_data_2)])

correlation_matrix <- cor(cor_data_2, method = "pearson",use = "pairwise.complete.obs")

# Plotting the heatmap
quartz()
heatmap(
  correlation_matrix,
  col = my_palette,
  symm = TRUE, # Display the heatmap symmetrically
  margins = c(10, 10), # Add margins around the heatmap
  main = "Watered", 
  Rowv = NA, 
  Colv = NA
)


#MONOCULTURES
cor_data = merged[1:432,]
cor_data = cor_data[cor_data[,"richness"] == 1,]
cor_data_2 = cbind(cor_data[6:13], cor_data[18:33])

numeric_cols <- sapply(cor_data_2, is.numeric)
cor_data_2[, !numeric_cols] <- lapply(cor_data_2[, !numeric_cols], as.numeric)
cor_data_2[,1:ncol(cor_data_2)] = scale(cor_data_2[,1:ncol(cor_data_2)])

correlation_matrix <- cor(cor_data_2, method = "pearson",use = "pairwise.complete.obs")

# Plotting the heatmap
quartz()
heatmap(
  correlation_matrix,
  col = my_palette,
  symm = TRUE, # Display the heatmap symmetrically
  margins = c(10, 10), # Add margins around the heatmap
  main = "Monocultures", 
  Rowv = NA, 
  Colv = NA
)

#POLYCULTURES
cor_data = merged[1:432,]
cor_data = cor_data[cor_data[,"richness"] == 12,]
cor_data_2 = cbind(cor_data[6:13], cor_data[18:33])

numeric_cols <- sapply(cor_data_2, is.numeric)
cor_data_2[, !numeric_cols] <- lapply(cor_data_2[, !numeric_cols], as.numeric)
cor_data_2[,1:ncol(cor_data_2)] = scale(cor_data_2[,1:ncol(cor_data_2)])

correlation_matrix <- cor(cor_data_2, method = "pearson",use = "pairwise.complete.obs")

# Plotting the heatmap
quartz()
heatmap(
  correlation_matrix,
  col = my_palette,
  symm = TRUE, # Display the heatmap symmetrically
  margins = c(10, 10), # Add margins around the heatmap
  main = "Polycultures", 
  Rowv = NA, 
  Colv = NA
)







