#Damla Cinoglu and Sarah Ortiz et al., 2025
#Checked by Damla Cinoglu in February 2025
####################################################
#SETUP: PACKAGES, DATA, DIRECTORY
####################################################

library(readxl)
library(plyr)
library(readr)
library(data.table)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggthemes)
library(ggpubr)
library(gridExtra)
library(devtools)
library(plotrix)
library(lme4)
library(gplots) 
library(ggrepel)
library(stats)
library(GGally)
library(ggcorrplot)
library(patchwork)
library(cowplot)
library(moments)
library(nlme)

setwd("/Users/damlacinoglu/Desktop/final_bfl_data")
merged = read.table(file = "MERGED_231122")

####################################################
#EDIT PCA DATA
####################################################

#add leaf C:N
merged = cbind(merged, c_N = merged$percC/ merged$percN)

#add leaf N by area
#convert LMA from kg/m2 to mg/m2 
merged = cbind(merged, percC_area = rep(NA, nrow(merged)), percN_area = rep(NA, nrow(merged)))
#LMA in kg/m2 to mg/m2 
#merged$LMA * 1000000, then divide by 1000000 again to convert to kg/m2

merged$percC = merged$percC/100
merged$percN = merged$percN/100

#in kg/m2
merged$percC_area = merged$percC * merged$LMA 
merged$percN_area = merged$percN * merged$LMA

#only include data from 2022
merged = merged[1:432,]

#Save the data to send to Sarah
write.csv(merged, "MERGED_231122_250326")



#removed COTI
merged = merged[merged[,"species"] != "COTI",]

#only focus on relevant traits 
merged = merged[,c(1,2,16,4,5,6,7,9,10,12,21,25,20,32,34,36)]
merged = na.omit(merged)

#TO MAKE THE CONSERVATIVE STUFF ACTUALLY CONSERVATIVE:
merged[,"LMA"] = -1 *merged[,"LMA"]
merged[,"d13C"]= -1 *merged[,"d13C"]
merged[,"c_N"]= -1 *merged[,"c_N"]
merged[,"root.tissue.density"]= -1 *merged[,"root.tissue.density"]
#merged[,"Depth.mm"]= -1 *merged[,"Depth.mm"]
merged[,"depth.cm"]= -1 *merged[,"depth.cm"]
merged[,"Average.Diameter.mm"]= -1 *merged[,"Average.Diameter.mm"]

####################################################
#FIGURE S1: PANEL A AND B
####################################################
diagnostics = F
plot_list = list()

for(i in c("growth", "species")) {
traits_data = merged

#get rid of NA entries
traits_data_clean = na.omit(traits_data)

traits_data_clean_clean = traits_data_clean[,-c(1,3)]

if(diagnostics == T) {
	# Pairwise comparisons to check linearity
scaled_data = traits_data_clean_clean %>%
  dplyr::select(-treatment) %>%
  dplyr::select(-richness) %>%
  dplyr::select(-species) %>%
  scale(center = TRUE, scale = TRUE) %>%
  as.data.frame()
# pairwise_plots = ggpairs(
#   scaled_data,
#   upper = list(continuous = wrap("smooth_lm", se = FALSE)),  # Regression line without confidence interval
#   lower = list(continuous = wrap("smooth_lm", se = FALSE)),  # Regression line without confidence interval
#   diag = list(continuous = "densityDiag")  # Density plot on the diagonal
# )

#check for multicollinearity?
#correlation_matrix = cor(traits_data_clean_clean %>% select(-treatment))
correlation_matrix = cor(scaled_data)
print(correlation_matrix)
ggcorrplot(correlation_matrix, hc.order = TRUE, type = "lower", lab = TRUE)
}

traits_data_clean_clean = traits_data_clean_clean[ , !(names(traits_data_clean_clean) %in% "percN")]

traits_data_clean_clean = traits_data_clean_clean[complete.cases(traits_data_clean_clean) & 
                                                     !apply(traits_data_clean_clean, 1, function(x) any(is.infinite(x))), ]
#the where is numeric is getting rid of the "treatment" column
traits_data_clean_clean$veg_height = as.numeric(traits_data_clean_clean[,"veg_height"])
pca_result = prcomp(traits_data_clean_clean[,4:ncol(traits_data_clean_clean)] , scale. = TRUE, center = TRUE)

#how much variation is being explained by two axes:
#pca_summary = summary(pca_result)
#explained_variance = pca_result$sdev^2 / sum(pca_result$sdev^2)
#variance_PC1 = explained_variance[1] * 100  # For the 1st principal component
#variance_PC2 = explained_variance[2] * 100  # For the 2nd principal component
#cat("Variance explained by PC1:", variance_PC1, "%\n")
#cat("Variance explained by PC2:", variance_PC2, "%\n")

# Extract the principal components
pca_data = as.data.frame(pca_result$x)

# Add species column to the PCA data
if(i == "species") {pca_data$species = traits_data_clean_clean$species}
if(i == "growth") {pca_data$species = traits_data_clean$growth.form}
pca_data$treatment = traits_data_clean_clean$treatment
pca_data$richness = traits_data_clean_clean$richness

# Extract standard deviations and rotation matrix
rotation_matrix = pca_result$rotation

std_dev = pca_result$sdev
label_offset = 3  # Adjust this value as needed

# Create a data frame for arrows
arrow_data = data.frame(
  PC1 = rep(0, ncol(rotation_matrix)),
  PC2 = rep(0, ncol(rotation_matrix)),
  arrow_x = std_dev[1] * rotation_matrix[, 1],
  arrow_y = std_dev[2] * rotation_matrix[, 2],
  column_name = rownames(rotation_matrix)
)

# Plot PCA with arrows and colored points
arrow_data[,"arrow_x"] = label_offset * arrow_data[,"arrow_x"]
arrow_data[,"arrow_y"] = label_offset * arrow_data[,"arrow_y"] 


arrow_data$column_name = c("LA","LMA","d13C","Height","Diameter","SRL","Depth","RTD","LC:N","LN/area")

pca_plot = ggplot(pca_data, aes(x = PC1, y = PC2, color = species)) +
  geom_point() +
  geom_segment(data = arrow_data,
               aes(x = PC1, y = PC2, xend = PC1 + arrow_x, yend = PC2 + arrow_y),
               arrow = arrow(length = unit(0.1, "inches")),
               color = "black") +
  geom_text(data = arrow_data,
            aes(label = column_name,
                x = PC1 + arrow_x + sign(arrow_x) * 0.2,
                y = PC2 + arrow_y + sign(arrow_y) * 0.2),
            hjust = -0.2, vjust = 0.5, size = 3, color = "black") +
            ylim(c(-2.5, 3.1))+xlim(c(-4,4))+
  theme_bw() +
  if(i == "species") {labs(x = "PC1", y = "PC2", color = "Species")}else {labs(x = "PC1", y = "PC2", color = "Functional Group")} 
  
  plot_list[[i]] = pca_plot
}

####################################################
#FIGURE 3: PANEL A / FIGURE S1: PANEL C
####################################################
  
  # Define the trait organ categorization
arrow_data$trait_organ = ifelse(arrow_data$column_name %in% c("LA", "LMA","LN/area","LC:N", "d13C", "Height"), 
                                 "Aboveground", "Belowground")
#arrow_data[,"arrow_x"] = label_offset * arrow_data[,"arrow_x"]
#arrow_data[,"arrow_y"] = label_offset * arrow_data[,"arrow_y"] 
 
# Update ggplot code
above_below = ggplot(pca_data, aes(x = PC1, y = PC2)) +
  geom_point(color = "grey") +  # All points in grey
  geom_segment(data = arrow_data,
               aes(x = 0, y = 0, 
                   xend = arrow_x, yend = arrow_y, 
                   color = trait_organ),
               arrow = arrow(length = unit(0.1, "inches")),
               size = 0.8) +  # Arrow size
  geom_text(data = arrow_data,
            aes(label = column_name,
                x = arrow_x + sign(arrow_x) * 0.2,
                y = arrow_y + sign(arrow_y) * 0.2),
            hjust = -0.2, vjust = 0.5, size = 3, color = "black") +  # Remove color from legend
  scale_color_manual(values = c("Aboveground" = "black", "Belowground" = "red")) +  # Custom colors for traits
  ylim(c(-2.5, 3.1)) + 
  xlim(c(-4, 4)) +
  labs(x = "PC1", y = "PC2", color = "Trait Organ") +  # Legend title for arrows
  theme_bw()

final_plot = plot_list$species + plot_list$growth + above_below

ggsave("figs1.png", final_plot, width = 33, height = 11, units = "cm", dpi = 300)

####################################################
#FIGURE 3: PANEL B AND C
####################################################
n_unique_points <- nrow(pca_data[
  pca_data[,"PC2"] > 3.1 | 
  pca_data[,"PC2"] < -2.5 | 
  pca_data[,"PC1"] < -4 | 
  pca_data[,"PC1"] > 4, 
])

colnames(pca_data)[11] = "Species"

# WATERING
watering_centroids = pca_data %>%
  group_by(Species, treatment) %>%
  summarize(across(c(PC1, PC2), mean, na.rm = TRUE), .groups = "drop")

watering_centroids[watering_centroids[,"treatment"] == "C","treatment"] = "N"
watering_centroids$treatment <- factor(watering_centroids$treatment, 
                                       levels = c("N", "W"))

water_arrow_data <- watering_centroids %>%
  pivot_wider(names_from = treatment, values_from = c(PC1, PC2)) %>%
  rename(x_start = PC1_N, y_start = PC2_N, x_end = PC1_W, y_end = PC2_W)
  
watering = ggplot(watering_centroids, aes(x = PC1, y = PC2, shape = treatment)) +
  # Grey PCA arrows (Plotted first to appear behind)
 geom_segment(data = arrow_data, inherit.aes = FALSE,
               aes(x = 0, y = 0, xend = arrow_x, yend = arrow_y), 
               color = "grey", 
               arrow = arrow(length = unit(0.1, "inches")),
               size = 0.8) +  
  # PCA arrow labels
    geom_text(data = arrow_data, inherit.aes = FALSE,
            aes(label = column_name,
                x = arrow_x + sign(arrow_x) * 0.2,
                y = arrow_y + sign(arrow_y) * 0.2),
            hjust = -0.2, vjust = 0.5, size = 3, color = "grey") +  
              geom_point(color = "black", size = 3) +
  # Species-colored arrows connecting circle-to-triangle
  geom_segment(data = water_arrow_data, inherit.aes = FALSE,
               aes(x = x_start, y = y_start, 
                   xend = x_end, yend = y_end, 
                   color = Species),  
               arrow = arrow(length = unit(0.1, "inches")),
               size = 0.8) +
  ylim(c(-2.5, 3.1)) + xlim(c(-4, 4)) +
  labs(x = "PC1", y = "PC2", shape = "Treatment") +
  theme_bw() +
  guides(shape = guide_legend(order = 1), 
         color = guide_legend(order = 2))

# RICHNESS
richness_centroids = pca_data %>%
  group_by(Species, richness) %>%
  summarize(across(c(PC1, PC2), mean, na.rm = TRUE), .groups = "drop")
richness_centroids$richness <- as.character(richness_centroids$richness)
richness_centroids[richness_centroids[,"richness"] == "1", "richness"] <- "M"
richness_centroids[richness_centroids[,"richness"] == "12", "richness"] <- "P"
richness_centroids$richness <- as.factor(richness_centroids$richness)

colnames(richness_centroids)[2] = "Treatment"

# Create arrow data for richness
richness_arrow_data <- richness_centroids %>%
  pivot_wider(names_from = richness, values_from = c(PC1, PC2)) %>%
  rename(x_start = PC1_M, y_start = PC2_M, 
         x_end = PC1_P, y_end = PC2_P)

# Richness Plot
richness = ggplot(richness_centroids, aes(x = PC1, y = PC2, shape = Treatment)) +
  # Grey PCA arrows (Plotted first to appear behind)
  geom_segment(data = arrow_data, inherit.aes = FALSE,
               aes(x = 0, y = 0, xend = arrow_x, yend = arrow_y), 
               color = "grey", 
               arrow = arrow(length = unit(0.1, "inches")),
               size = 0.8) +
  # PCA arrow labels
  geom_text(data = arrow_data, inherit.aes = FALSE,
            aes(label = column_name,
                x = arrow_x + sign(arrow_x) * 0.2,
                y = arrow_y + sign(arrow_y) * 0.2),
            hjust = -0.2, vjust = 0.5, size = 3, color = "grey") +
  # Black points for species centroids
  geom_point(color = "black", size = 3) +
  # Species-colored arrows connecting circle-to-triangle
  geom_segment(data = richness_arrow_data, inherit.aes = FALSE,
               aes(x = x_start, y = y_start, 
                   xend = x_end, yend = y_end, 
                   color = Species),  
               arrow = arrow(length = unit(0.1, "inches")),
               size = 0.8) +
  ylim(c(-2.5, 3.1)) + xlim(c(-4, 4)) +
  labs(x = "PC1", y = "PC2", shape = "Treatment") +
  theme_bw() +
  guides(shape = guide_legend(order = 1), 
         color = guide_legend(order = 2))

final_plot = above_below + watering + richness  

ggsave("fig3.png", final_plot, width = 33, height = 11, units = "cm", dpi = 300)
 
####################################################
#FIGURE S2: LEAF + ROOT PCA FIGURE FOR SUPPLEMENT 
####################################################
traits_data = merged

traits_data_clean = na.omit(traits_data)

traits_data_clean_clean = traits_data_clean[,-c(1,3)]

traits_data_clean_clean = traits_data_clean_clean[ , !(names(traits_data_clean_clean) %in% "percN")]

######################Leaf and Root trait PCAs separately:
pca_plots = list()

for(i in c("leaf", "root")) {
	if(i == "leaf") {traits_data_clean_clean = traits_data_clean_clean[,c("species","richness","treatment","LA","LMA","d13C","veg_height" ,"c_N","percN_area")]
		traits_data_clean_clean$veg_height = as.numeric(traits_data_clean_clean$veg_height )}
 	if(i == "root") {
 		traits_data_clean_clean = traits_data_clean[,-c(1,3)]

traits_data_clean_clean = traits_data_clean_clean[ , !(names(traits_data_clean_clean) %in% "percN")]

 		traits_data_clean_clean = traits_data_clean_clean[,c("species","richness","treatment","Average.Diameter.mm","srl","Depth.mm","root.tissue.density")]}
	
traits_data_clean_clean = traits_data_clean_clean[complete.cases(traits_data_clean_clean) & 
                                                     !apply(traits_data_clean_clean, 1, function(x) any(is.infinite(x))), ]
#the where is numeric is getting rid of the "treatment" column
pca_result = prcomp(traits_data_clean_clean[,4:ncol(traits_data_clean_clean)] , scale. = TRUE, center = TRUE)

#how much variation is being explained by two axes:
pca_summary = summary(pca_result)
print(pca_summary)
explained_variance = pca_result$sdev^2 / sum(pca_result$sdev^2)
print(explained_variance)
variance_PC1 = explained_variance[1] * 100  # For the 1st principal component
variance_PC2 = explained_variance[2] * 100  # For the 2nd principal component
cat("Variance explained by PC1:", variance_PC1, "%\n")
cat("Variance explained by PC2:", variance_PC2, "%\n")

# Extract the principal components
pca_data = as.data.frame(pca_result$x)

# Add species column to the PCA data
pca_data$species = traits_data_clean_clean$species
pca_data$treatment = traits_data_clean_clean$treatment
pca_data$richness = traits_data_clean_clean$richness

# Extract standard deviations and rotation matrix
std_dev = pca_result$sdev
rotation_matrix = pca_result$rotation

# Create a data frame for arrows
arrow_data = data.frame(
  PC1 = rep(0, ncol(rotation_matrix)),
  PC2 = rep(0, ncol(rotation_matrix)),
  arrow_x = std_dev[1] * rotation_matrix[, 1],
  arrow_y = std_dev[2] * rotation_matrix[, 2],
  column_name = rownames(rotation_matrix))

scale_factor = 3
arrow_data$arrow_x = arrow_data$arrow_x * scale_factor
arrow_data$arrow_y = arrow_data$arrow_y * scale_factor

 if(i == "leaf") {arrow_data[,"column_name"] = c("LA","LMA","d13C","Height","LC:N","LN/area")}
 if(i == "root") {arrow_data[,"column_name"] = c("Diameter","SRL","Depth","RTD")}
 
pca_plot = ggplot(pca_data, aes(x = PC1, y = PC2, color = species)) +
  geom_point() +
  geom_segment(data = arrow_data,
               aes(x = PC1, y = PC2, xend = PC1 + arrow_x, yend = PC2 + arrow_y),
               arrow = arrow(length = unit(0.1, "inches")),
               color = "black") +
  geom_text(data = arrow_data,
            aes(label = column_name,
                x = PC1 + arrow_x + sign(arrow_x) * 0.2,
                y = PC2 + arrow_y + sign(arrow_y) * 0.2),
            hjust = -0.2, vjust = 0.5, size = 3, color = "black") +
  #labs(x = "PC1", y = "PC2", color = "Functional Group") +
  labs(x = "PC1", y = "PC2", color = "Species") +
  ylim(c(-2.8, 3.1)) + xlim(c(-4, 4)) +
  theme_bw()
  
 if(i == "leaf") {leaf_pca = pca_data}
 if(i == "root") {root_pca = pca_data}
 
 pca_plots[[i]] = pca_plot
}

final_plot = pca_plots$leaf + pca_plots$root

# Fit the linear model
model = lm(root_pca[root_pca [,"species"] %in% c("LUTE","MOCI","GAPU","CEAM","PHCO","MOPU","IPRU","DEIL"),"PC1"] ~ leaf_pca[leaf_pca [,"species"] %in% c("LUTE","MOCI","GAPU","CEAM","PHCO","MOPU","IPRU","DEIL"),"PC1"])
summary(model)

model = lm(root_pca[root_pca [,"species"] %in% c("ARPU","BOCU","BOGR"),"PC1"] ~ leaf_pca[leaf_pca [,"species"] %in% c("ARPU","BOCU","BOGR"),"PC1"])
summary(model)

# Define species groups
non_grass_species = c("LUTE", "MOCI", "GAPU", "CEAM", "PHCO", "MOPU", "IPRU", "DEIL")
grass_species = c("ARPU", "BOCU", "BOGR")

# Filter data
leaf_non_grass = leaf_pca %>% filter(species %in% non_grass_species)
root_non_grass = root_pca %>% filter(species %in% non_grass_species)

leaf_grass = leaf_pca %>% filter(species %in% grass_species)
root_grass = root_pca %>% filter(species %in% grass_species)

# Combine into data frames
pca_non_grass = data.frame(
  species = leaf_non_grass$species,
  Leaf_PC1 = leaf_non_grass$PC1,
  Root_PC1 = root_non_grass$PC1)

pca_grass = data.frame(
  species = leaf_grass$species,
  Leaf_PC1 = leaf_grass$PC1,
  Root_PC1 = root_grass$PC1)

# Plot for non-grasses
plot1 = ggplot(pca_non_grass, aes(x = Leaf_PC1, y = Root_PC1)) +
  geom_point(color = "black", size = 2, alpha = 0.8) +
  geom_smooth(method = "lm", color = "red", se = FALSE) +
  labs(
    x = "Leaf PC1",
    y = "Root PC1",
    title = "Root PC1 vs. Leaf PC1 (Non-Grasses)"
  ) +theme_minimal() +
  theme_bw()

# Plot for grasses
plot2 = ggplot(pca_grass, aes(x = Leaf_PC1, y = Root_PC1)) +
  geom_point(color = "black", size = 2, alpha = 0.8) +
  #geom_smooth(method = "lm", color = "red", se = FALSE) +
  labs(
    x = "Leaf PC1",
    y = "Root PC1",
    title = "Root PC1 vs. Leaf PC1 (Grasses)"
  ) +theme_minimal() +
  theme_bw()

top_row <- plot_grid(final_plot, ncol = 1)  
bottom_row <- plot_grid(plot2, plot1, ncol = 2)  
final_plot = plot_grid(top_row, bottom_row, ncol = 1)  

ggsave("figs2.png", final_plot, width = 23, height = 23, units = "cm", dpi = 300)

####################################################
#FIGURE S3: ABIOTIC VARIABLE PLOTS
####################################################
setwd("/Users/damlacinoglu/Desktop/final_bfl_data")
abiotic_merged = read.table( file = "ABIOTIC_240401*.txt")
abiotic_merged = abiotic_merged %>% distinct(plot, .keep_all = TRUE)

nrow(abiotic_merged)
unique(length(abiotic_merged$temp))
unique(length(abiotic_merged$moist))
unique(length(abiotic_merged$nh4))
unique(length(abiotic_merged$light_atten))


temp_summary = abiotic_merged %>%
  drop_na(temp) %>% 
  group_by(treatment, richness) %>%
  summarise_at(vars(temp), list(sd = sd,
                                 se = std.error,  
                                 temp = mean))
                                 
moist_summary = abiotic_merged %>%
  drop_na(moist) %>% 
  group_by(treatment, richness) %>%
  summarise_at(vars(moist), list(sd = sd,
                                 se = std.error,  
                                 temp = mean))

nh4_summary = abiotic_merged %>%
  drop_na(nh4) %>%
  #mutate(log_nh4 = log(nh4)) %>% # Calculate the natural log of nh4
  group_by(treatment, richness) %>%
  #summarise_at(vars(nh4, log_nh4), list(sd = sd, 
   summarise_at(vars(nh4), list(sd = sd, 
  										se = std.error, 
  										temp = mean))
  										
 abiotic_merged[,"light_atten"] = 
  abiotic_merged[,"light_canopy"] - ((abiotic_merged[,"light_ground1"]+abiotic_merged[,"light_ground2"])/2) 

light_summary = abiotic_merged %>%
  drop_na(light_atten) %>%
  #mutate(log_nh4 = log(nh4)) %>% # Calculate the natural log of nh4
  group_by(treatment, richness) %>%
  #summarise_at(vars(nh4, log_nh4), list(sd = sd, 
   summarise_at(vars(light_atten), list(sd = sd, 
  										se = std.error, 
  										temp = mean))
  										
# Function to create a plot for a specific variable
create_plot = function(summary_data, y_label, y_limits = c(0,1), reverse_y = FALSE) {
  ymax = max(summary_data$temp + summary_data$se) + 3
  
  p <- ggplot(summary_data, aes(x = factor(richness), y = temp, color = treatment)) +
    geom_point(size = 4, position = position_dodge(width = 0.25)) +
    geom_errorbar(aes(ymin = temp - se, ymax = temp + se), width = 0.5, position = position_dodge(width = 0.25)) +
    ylab(y_label) +  
    xlab('') + 
    theme(text = element_text(size = 12)) + 
    theme_bw(base_size = 12) + 
    theme(
      strip.text = element_text(size = 10),
      panel.background = element_rect(fill = 'transparent'),
      plot.background = element_rect(fill = 'transparent', color = NA),
      legend.background = element_rect(fill = 'transparent'),
      legend.box.background = element_rect(fill = 'transparent'),
      legend.position = 'bottom'
    ) +
    scale_color_manual(values= c('skyblue2', 'skyblue4'), labels=c('N', 'W')) +
    labs(color = " ") +
    scale_x_discrete(labels = c("M", "P"))
  
  # Apply y-axis transformation
  if (reverse_y) {
    p <- p + scale_y_reverse(limits = rev(y_limits))
  } else {
    p <- p + ylim(y_limits)
  }
  
  return(p)
}

# Create individual plots
plot_temp = create_plot(temp_summary, "Soil temperature (C)", c(50,59))
plot_moist = create_plot(moist_summary, "Soil moisture (%)",c(0.08,0.17))
plot_nh4 = create_plot(nh4_summary, "Soil NH4 (ug/g)",c(0.5,4.9))
plot_light = create_plot(light_summary, "Light attenuation (mol m^-2 d^-1)", c(300,1500), reverse_y = TRUE)

# Combine plots using cowplot
combined_plot = plot_grid(
  plot_temp + theme(legend.position = "none"), 
  plot_moist + theme(legend.position = "none"), 
  plot_nh4 + theme(legend.position = "none"), 
    plot_light + theme(legend.position = "none"), 
  nrow = 1)

ggsave("figs3_new.png", combined_plot, width = 24, height = 9, units = "cm", dpi = 300)

abiotic_merged = cbind(abiotic_merged, block = rep(NA, nrow(abiotic_merged)))
for(i in 1:nrow(abiotic_merged)){ 
  plot = abiotic_merged[i,"plot"] 
if(grepl("A",plot)) {abiotic_merged[i,"block"] = "A"}
  if(grepl("B",plot)) {abiotic_merged[i,"block"] = "B"}
  if(grepl("C",plot)) {abiotic_merged[i,"block"] = "C"}
}

kurtosis(abiotic_merged$temp, na.rm = TRUE)
skewness(abiotic_merged$temp, na.rm = TRUE)
srl.wat.aov = aov(temp ~as.factor(richness)*treatment + Error(block), data = abiotic_merged)
summary(srl.wat.aov) 

kurtosis(abiotic_merged$moist, na.rm = TRUE)
skewness(abiotic_merged$moist, na.rm = TRUE)
srl.wat.aov = aov(moist ~as.factor(richness)*treatment + Error(block), data = abiotic_merged)
summary(srl.wat.aov) 

kurtosis(abiotic_merged$nh4, na.rm = TRUE)
skewness(abiotic_merged$nh4, na.rm = TRUE)
abiotic_merged$nh4 = log(abiotic_merged$nh4+1)
srl.wat.aov = aov(nh4 ~ as.factor(richness)*treatment + Error(block), data = abiotic_merged)
summary(srl.wat.aov)

kurtosis(abiotic_merged$light_atten, na.rm = TRUE)
skewness(abiotic_merged$light_atten, na.rm = TRUE)
srl.wat.aov = aov(light_atten ~richness*treatment + Error(block), data = abiotic_merged)
summary(srl.wat.aov)

####################################################
#FIGURE S3: COMMUNITY BIOMASS PLOT
####################################################
setwd("/Users/damlacinoglu/Desktop/sarah trait paper check")
biomass = data.frame(read_excel("BIOMASS_2022.xlsx"))
biomass = biomass[biomass[,"Quadrat"] == "Open",]

new_biomass = data.frame(unique(biomass$Plot.name))
new_biomass <- na.omit(new_biomass)
new_biomass = cbind(new_biomass, total_biomass = rep(NA, nrow(new_biomass)))
for(i in 1:nrow(new_biomass)) {
		plot = new_biomass[i,"unique.biomass.Plot.name."]
new_biomass[i,"total_biomass"] = sum(na.omit(as.numeric(biomass[biomass[,"Plot.name"] == plot & biomass[,"Quadrat"] == "Open","Total.weight..g."])))
}
new_biomass = cbind(new_biomass, richness = rep(NA, nrow(new_biomass)))
new_biomass = cbind(new_biomass, watering = rep(NA, nrow(new_biomass)))
new_biomass = cbind(new_biomass, block = rep(NA, nrow(new_biomass)))
for(i in 1:nrow(new_biomass)){ 
  plot = new_biomass[i,"unique.biomass.Plot.name."] 
if(grepl("A",plot)) {new_biomass[i,"block"] = "A"}
  if(grepl("B",plot)) {new_biomass[i,"block"] = "B"}
  if(grepl("C",plot)) {new_biomass[i,"block"] = "C"}
}
new_biomass = new_biomass %>%
  mutate(watering = ifelse(grepl("W", unique.biomass.Plot.name.), "W", "C"))
  for(i in 1:nrow(new_biomass)){ 
  plot = new_biomass[i,"unique.biomass.Plot.name."] 
if(grepl("MON",plot)) {new_biomass[i,"richness"] = 1} else {
	if(grepl("2",plot)) {new_biomass[i,"richness"] = 2}
  if(grepl("4",plot)) {new_biomass[i,"richness"] = 4}
	  if(grepl("6",plot)) {new_biomass[i,"richness"] = 6}
	  if(grepl("12",plot)) {new_biomass[i,"richness"] = 12}
}
}

#Save biomass data to send to sarah:
write.csv(new_biomass, "biomass_summary_2022.csv", row.names = FALSE)

# Filter and modify richness values for monoculture and polyculture
new_biomass_just = new_biomass %>%
  filter(richness %in% c(1, 4, 6, 12)) 

# Summarize the data
biomass_summary = new_biomass_just %>%
  drop_na(total_biomass) %>% 
  group_by(watering, richness) %>%
  summarise(
    sd = sd(total_biomass),
    se = std.error(total_biomass),
    temp = mean(total_biomass),
    .groups = "drop"
  )
  
# Modify the create_plot() function to use richness_grouped
create_plot2 = function(summary_data, y_label, y_limits = c(0, 1)) {
  ymax = max(summary_data$temp + summary_data$se, na.rm = TRUE) + 3

  ggplot(summary_data, aes(x = factor(richness), y = temp, color = watering)) +
    geom_point(size = 4, position = position_dodge(width = 0.25)) +
    geom_errorbar(aes(ymin = temp - se, ymax = temp + se), 
                  width = 0.5, position = position_dodge(width = 0.25)) +
    ylab(y_label) +  
    xlab('') + 
    ylim(y_limits) + 
    theme(text = element_text(size = 12)) + 
    theme_bw(base_size = 12) + 
    theme(
      strip.text = element_text(size = 10),
      panel.background = element_rect(fill = 'transparent'),
      plot.background = element_rect(fill = 'transparent', color = NA),
      legend.background = element_rect(fill = 'transparent'),
      legend.box.background = element_rect(fill = 'transparent'),
      legend.position = 'bottom'
    ) +
    scale_color_manual(values = c('skyblue2', 'skyblue4'), labels = c('N', 'W')) +
    labs(color = " ") +
    scale_x_discrete(labels = c("1" = "M", "4" = "P = 4", "6" = "P = 6", "12" = "P = 12"))
}
  
plot_biomass = create_plot2(biomass_summary, "Community aboveground biomass (g)", c(25,250))

kurtosis(log(new_biomass_just $total_biomass+1), na.rm = TRUE)
skewness(log(new_biomass_just $total_biomass+1), na.rm = TRUE)
  srl.wat.aov = aov(log(new_biomass_just $total_biomass+1) ~ as.factor(richness)*watering + Error(block), data = new_biomass_just)
summary(srl.wat.aov) 
    srl.wat.aov = aov(log(new_biomass_just $total_biomass+1) ~ richness*watering + Error(block), data = new_biomass_just)
summary(srl.wat.aov) 

library(emmeans)
emmeans_results <- emmeans(srl.wat.aov, pairwise ~ richness * watering, adjust = "tukey")
summary(emmeans_results)

ggsave("figs3.png", plot_biomass, width = 10, height = 10, units = "cm", dpi = 300)

####################################################
#MANUAL ADDITIONS TO PLOTS
####################################################
#Nonwatered + Watered legend at the botom for Figure S3
#Panel A, B, C, D markers at the top left corner of plots
#Add carets to traits in Figure 3, Figure S1, Figure S2?
#Figure S3 all panels and Figure S4, all have significance and model values noted on top right corner

####################################################
#ALTERNATIVE FIGURE S3: NEW ABIOTIC VARIABLE PLOTS
####################################################
setwd("/Users/damlacinoglu/Desktop/sarah trait paper check")
abiotic_merged = read.table( file = "ABIOTIC_250321")
abiotic_merged = abiotic_merged %>% distinct(plot, .keep_all = TRUE)

temp_merged <- abiotic_merged %>%
  select(individual, species, plot, richness, treatment, contains("temp")) %>%
  pivot_longer(cols = contains("temp"),
               names_to = "temp_variable",
               values_to = "temp_value") %>%
  mutate(
    season = sub("temp_(mid|end)_.*", "\\1", temp_variable),
    year = sub(".*_(\\d{2})$", "\\1", temp_variable)
  ) %>%
  select(individual, species, plot, richness, treatment, year, season, temp_value)
temp_summary = temp_merged %>%
  drop_na(temp_value) %>% 
  group_by(treatment, richness) %>%
  summarise_at(vars(temp_value), list(sd = sd,
                                 se = std.error,  
                                 temp = mean))

moist_merged <- abiotic_merged %>%
  select(individual, species, plot, richness, treatment, contains("moist")) %>%
  pivot_longer(cols = contains("moist"),
               names_to = "moist_variable",
               values_to = "moist_value") %>%
  mutate(
    season = sub("moist_(mid|end)_.*", "\\1", moist_variable),
    year = sub(".*_(\\d{2})$", "\\1", moist_variable)
  ) %>%
  select(individual, species, plot, richness, treatment, year, season, moist_value)                       
moist_summary = moist_merged %>%
  drop_na(moist_value) %>% 
  group_by(treatment, richness) %>%
  summarise_at(vars(moist_value), list(sd = sd,
                                 se = std.error,  
                                 temp = mean))

nh4_summary = abiotic_merged %>%
  drop_na(nh4) %>%
  #mutate(log_nh4 = log(nh4)) %>% # Calculate the natural log of nh4
  group_by(treatment, richness) %>%
  #summarise_at(vars(nh4, log_nh4), list(sd = sd, 
   summarise_at(vars(nh4), list(sd = sd, 
  										se = std.error, 
  										temp = mean))
  										
 abiotic_merged[,"light_atten_23"] = 
  ((abiotic_merged[,"light_ground1_end_23"]+abiotic_merged[,"light_ground2_end_23"])/2) - abiotic_merged[,"light_canopy_end_23"]
 abiotic_merged[,"light_atten_24"] = 
  ((abiotic_merged[,"light_ground1_end_24"]+abiotic_merged[,"light_ground2_end_24"])/2) - abiotic_merged[,"light_canopy_end_24"]
  light_merged <- abiotic_merged %>%
  select(individual, species, plot, richness, treatment, 
         light_atten_23, light_atten_24) %>%
  pivot_longer(cols = c(light_atten_23, light_atten_24),
               names_to = "light_variable",
               values_to = "light_value") %>%
  mutate(year = sub(".*_(\\d{2})$", "\\1", light_variable)) %>%
  select(individual, species, plot, richness, treatment, year, light_value)
light_summary = light_merged %>%
  drop_na(light_value) %>%
  #mutate(log_nh4 = log(nh4)) %>% # Calculate the natural log of nh4
  group_by(treatment, richness) %>%
  #summarise_at(vars(nh4, log_nh4), list(sd = sd, 
   summarise_at(vars(light_value), list(sd = sd, 
  										se = std.error, 
  										temp = mean))
  										
# Create individual plots
plot_temp = create_plot(temp_summary, "Soil temperature (C)", c(31,43))
plot_moist = create_plot(moist_summary, "Soil moisture (%)",c(0.08,0.17))
plot_nh4 = create_plot(nh4_summary, "Soil NH4 (ug/g)",c(0.5,4.9))
plot_light = create_plot(light_summary, "Light attenuation (mol m^-2 d^-1)",c(-1400, -400))

# Combine plots using cowplot
combined_plot = plot_grid(
  plot_temp + theme(legend.position = "none"), 
  plot_moist + theme(legend.position = "none"), 
  plot_nh4 + theme(legend.position = "none"), 
    plot_light + theme(legend.position = "none"), 
  nrow = 1)

ggsave("figs3_new.png", combined_plot, width = 24, height = 9, units = "cm", dpi = 300)

temp_merged = cbind(temp_merged, block = rep(NA, nrow(temp_merged)))
for(i in 1:nrow(temp_merged)){ 
  plot = temp_merged[i,"plot"] 
if(grepl("A",plot)) {temp_merged[i,"block"] = "A"}
  if(grepl("B",plot)) {temp_merged[i,"block"] = "B"}
  if(grepl("C",plot)) {temp_merged[i,"block"] = "C"}
}
srl.wat.lmer = lmer(temp_value ~ richness * treatment + 
                      (1 | block) + (1 | year) + (1 | season), 
                      data = temp_merged)
summary(srl.wat.lmer)

moist_merged = cbind(moist_merged, block = rep(NA, nrow(moist_merged)))
for(i in 1:nrow(moist_merged)){ 
  plot = moist_merged[i,"plot"] 
if(grepl("A",plot)) {moist_merged[i,"block"] = "A"}
  if(grepl("B",plot)) {moist_merged[i,"block"] = "B"}
  if(grepl("C",plot)) {moist_merged[i,"block"] = "C"}
}
srl.wat.lmer = lmer(moist_value ~ richness * treatment + 
                      (1 | block) + (1 | year) + (1 | season)+ (1 | plot), 
                      data = moist_merged)
summary(srl.wat.lmer)

abiotic_merged$nh4 = log(abiotic_merged$nh4+1)
srl.wat.aov = aov(nh4 ~richness*treatment + Error(block), data = abiotic_merged)
summary(srl.wat.aov)

kurtosis(abiotic_merged$light_atten, na.rm = TRUE)
skewness(abiotic_merged$light_atten, na.rm = TRUE)
light_merged = cbind(light_merged, block = rep(NA, nrow(light_merged)))
for(i in 1:nrow(light_merged)){ 
  plot = light_merged[i,"plot"] 
if(grepl("A",plot)) {light_merged[i,"block"] = "A"}
  if(grepl("B",plot)) {light_merged[i,"block"] = "B"}
  if(grepl("C",plot)) {light_merged[i,"block"] = "C"}
}
srl.wat.lmer = lmer(light_value ~ richness * treatment + 
                      (1 | block) + (1 | year) , 
                      data = light_merged)
summary(srl.wat.lmer)




