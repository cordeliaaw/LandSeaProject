# # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # 
# # #                             # # # 
# # #                             # # # 
# # #      Land Sea Project       # # # 
# # #                             # # # 
# # #                             # # # 
# # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # 





# # # # #  Description  # # # # # ----
#
#  The Land Sea Project analyses sediment cores taken across the Gulf of Bothnia 
#  in spring 2024 by A. Rouillard (UMF, Sweden). 
#  Methods include VNIR spectroscopy, XRF elemental analysis, 210Pb dating, 
#  C/N stable isotope analysis. 
#  This project aims to look at phytoplankton (using Chla as an indicator) trends 
#  including diatoms (using Si/detritals ratio as an indicator) in relation to
#  browning (inferred water column TOC). 
#  We look at riverine influence (coastal vs offshore) and at the different basins 
#  (Bothnian bay vs Bothnian Sea). 
# 
# # # # # # # # # # # # # # # # # # # #





###
### Project supervised by Dr. A. Rouillard and Prof. N. Kamenos (UMF)
###
### Script by A. Weisenfeld (MSc intern, Sorbonne University)
###





# # # # #  Resources  # # # # # 
#
# https://tombishop1.github.io/itraxBook/using-ordr-for-pca.html
# Pedersen, 2019
# https://github.com/eric-pedersen/mixed-effect-gams  
#
# # # # # # # # # # # # # # # # 





# # # # #  Overview  # # # # # ----
#
#  1 - Preparation
#     1.1 - Loading packages
#     1.2 - Environment
#     1.3 - Tidying data
#     1.4 - Quality control
# 
#  2 - Spatial variability
#     2.1 - Data exploration
#     2.2 - Monitoring vs phytoplankton
#     2.3 - Phytoplankton distribution
#     2.4 - Water colour
#     2.5 - Redox conditions
#
#  3 - Temporal variability
#     3.1 - Phytoplankton
#     3.2 - Browning
#     3.3 - Redox conditions
#
#  4 - Modeling
#     4.1 - TOC vs total phytoplankton
#     4.2 - TOC vs diatoms
#
#  5 - Misc
#     5.1 - Age-depth profiles
#     5.2 - Lead time series
#     5.3 - Testing universal model (TOC)
#     5.4 - Diagenetic layer
#     5.5 - Elemental pie chart
#     5.6 - Example code for saving plots (as pdf)
# 
# # # # # # # # # # # # # # # # # # 





# 1 Preparation ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

### 1.1 Loading packages ----

# install needed packages using:
# install.packages("")

# load packages

library(tidyverse) # includes ggplot2
library(tidypaleo) # for creating stratigraphic diagrams
library(patchwork) # to create more complex plots
library(forcats) # for reordering factor levels
library(ordr) # for PCA
library(compositions) # PCA
library(viridis) # colour palette
library(ggnewscale) # adding 2 scales for PCA plots
library(mgcv) # time series modeling (gams)
library(gratia) # graphics for gams 
library(dplyr)
library(stringi) # characters
library(RColorBrewer) # colour palette
library(tidyr)
library(ggpattern) # for adding elements to ggplots
library(GGally)
library(car)   # for leveneTest
library(rstatix) # for pairwise tests and significance letters
library(ggpubr) # for easy annotation








### 1.2 Environment ----


rm(list=ls()) # clear memory of all variables

setwd("~/Documents/Sorbonne/Master/Stage M2/data") # set working directory

quartz() # open plot in separate rsession window (for mac)
# use x11() for windows



# import data

XRF <- read.csv("XRF_complete_data.csv")
View(XRF)
# XRF contains elemental data (proportion of elements measured) + references

VNIR <- read.csv("VNIR_complete_data.csv")
View(VNIR)
# VNIR contains Chla data and TOC data along with spectral data + references

Base <- read.csv("LandSease_Base.csv")
View(Base)
# Base contains additional characteristics (water depth, moisture, etc.)

UniModel <- read.csv("TOC_UniversalModel.csv")
View(UniModel)
# UniModel contains model used to infer water column TOC (spectral data), data from northern lakes
# Meyer-Jacob et al, 2017

monitoring <- read.csv("UMF_monitoring_data.csv")
# monitoring contains monitoring data for the pas 3 decades
# here we focus on total phytoplankton biomass, which is the closest to our chl-a indicator



### 1.3 Tidying data ----

    #### XRF ----
# add column to XRF dataset to order from North to South

XRF <- XRF %>%
  mutate(site_order = case_when(
    Site == "Kal0" ~ 1,  
    Site == "Kal1" ~ 2,
    Site == "Kal2" ~ 3,
    Site == "Ran0" ~ 4,
    Site == "Ran1" ~ 5,
    Site == "Ran2" ~ 6,
    Site == "SE17" ~ 7,
    Site == "A5" ~ 8,
    Site == "A13" ~ 9,
    Site == "SE1" ~ 10,
    Site == "Ore1" ~ 11,
    Site == "Ore2" ~ 12,
    Site == "GA1" ~ 13,
    Site == "Gav1" ~ 14,
    Site == "Gav2" ~ 15,
    Site == "C3" ~ 16,
    Site == "C14" ~ 17,
    Site == "SE3" ~ 18,
    Site == "Ljs0" ~ 19,
    Site == "Ljs1" ~ 20,
    Site == "Ljs2" ~ 21,
    Site == "Dal2" ~ 22, # Dal2 excluded later on because insufficient data
    Site == "SE4" ~ 23,
    TRUE ~ NA_integer_  # fallback in case of unknown site names
  ))


# move column to position 7
XRF <- XRF %>%
  relocate(site_order, .before = 7)


# add ratios, adapt as necessary

XRF <- XRF %>%
  mutate(Mn_Fe = Mn / Fe)

XRF <- XRF %>%
  mutate(Ca_Sr = Ca / Sr)

XRF <- XRF %>%
  mutate(Si_Al = SiO2 / Al2O3)

XRF <- XRF %>%
  mutate(Si_Ti = SiO2 / Ti)

XRF <- XRF %>%
  mutate(Ca_Fe = Ca / Fe)

XRF <- XRF %>%
  mutate(Ca_Ti = Ca / Ti)

XRF <- XRF %>%
  mutate(Mn_Ti = Mn / Ti)

XRF <- XRF %>%
  mutate(Mn_Al = Mn / Al2O3)

XRF <- XRF %>%
  mutate(Mn_Ca = Mn / Ca)

XRF <- XRF %>%
  mutate(S_Al = S / Al2O3)

XRF <- XRF %>%
  mutate(S_Fe = S / Fe)

XRF <- XRF %>%
  mutate(Zr_Rb = Zr / Rb)

XRF <- XRF %>%
  mutate(Zr_Ti = Zr / Ti)

XRF <- XRF %>%
  mutate(Si_detrital = SiO2 / (Al2O3 + Rb + Ti + Zr))

XRF <- XRF %>%
  mutate(P_Fe = P / Fe)

XRF <- XRF %>%
  mutate(P_Al = P / Al2O3)


# XRF is the raw data, compiling the reference materials (use for quality control (QC)) and sediment data
# XRF_avg keeps only the sediment data and averages the replicates

# averaging replicates
XRF_avg <- XRF %>%
  filter(Sample_type == "Sediment") %>% # filter to keep only relevant sediment samples (not references)
  group_by(Site, Interval_top) %>%  # groups by core and interval (i.e. replicate definition: when equal, is replicate)
  summarise(across(where(is.numeric), \(x) mean(x, na.rm = TRUE)), .groups = "drop") # averages all numeric columns within those groups, then ungroups the result


# XRF_ref keeps only the references for QC
XRF_ref <- XRF %>%
  filter(Sample_type == "Reference") %>%
  mutate(order = row_number())  # approximate sequence

XRF_ref <- XRF_ref %>%
  relocate(order, .before = 9)





    #### VNIR ----
# add column to VNIR dataset to order from North to South

VNIR <- VNIR %>%
  mutate(site_order = case_when(
    Site == "Kal0" ~ 1,  
    Site == "Kal1" ~ 2,
    Site == "Kal2" ~ 3,
    Site == "Ran0" ~ 4,
    Site == "Ran1" ~ 5,
    Site == "Ran2" ~ 6,
    Site == "SE17" ~ 7,
    Site == "A5" ~ 8,
    Site == "A13" ~ 9,
    Site == "SE1" ~ 10,
    Site == "Ore1" ~ 11,
    Site == "Ore2" ~ 12,
    Site == "GA1" ~ 13,
    Site == "Gav1" ~ 14,
    Site == "Gav2" ~ 15,
    Site == "C3" ~ 16,
    Site == "C14" ~ 17,
    Site == "SE3" ~ 18,
    Site == "Ljs0" ~ 19,
    Site == "Ljs1" ~ 20,
    Site == "Ljs2" ~ 21,
    Site == "Dal2" ~ 22, # Dal2 excluded later on because insufficient data
    Site == "SE4" ~ 23,
    TRUE ~ NA_integer_  # fallback in case of unknown site names
  ))


# move column to position 6
VNIR <- VNIR %>%
  relocate(site_order, .before = 6)


# data correction for chla
VNIR <- VNIR %>%
  mutate(Chla_corr = VRS_chla + 0.003574034) # move all values up so that there are no negative values left (ie. lowest value = 0)


# move column to position 10
VNIR <- VNIR %>%
  relocate(Chla_corr, .before = 10)


# VNIR is the raw data, compiling the reference materials (use for quality control (QC)) and sediment data
# VNIR_avg keeps only the sediment data and averages the replicates

# averaging replicates

VNIR_avg <- VNIR %>%
  filter(Sample_type == "Sediment") %>%
  group_by(Site, Depth) %>%  # groups by core and depth (i.e. replicate definition)
  summarise(across(where(is.numeric), \(x) mean(x, na.rm = TRUE)), .groups = "drop") 
# averages all numeric columns within those groups, then ungroups the result


# VNIR_ref keeps only the references for QC

VNIR_ref <- VNIR %>%
  filter(Sample_type == "Reference") %>%
  mutate(order = row_number())  # approximate sequence



    #### Base ----

# removing accents and other special characters from Station_code column 

Base <- Base %>%
  mutate(Station_code = Station_code %>%
           stri_trans_general("Latin-ASCII") %>%  # Remove accents
           str_remove_all("[-_]"))                # Remove - and _





    #### Merging data ----

# merge XRF and VNIR data first 

merged_data <- XRF_avg %>%
  select(-site_order) %>%
  inner_join(
    VNIR_avg %>% select(-site_order),
    by = c("Site", "Depth")
  )
# merged data doesn't have Dal2
# site_order is removed, recreated later
# if site_order columns are identical in both XRF and VNIR datasets then keep site_order column when merging


# add Base dataset
merged_data <- merged_data %>%
  left_join(
    Base %>% select(-Site) %>%
      filter(Core == "C1"),
    by = c("Site" = "Station_code", "Depth" = "Mean_depth") # same columns have different names in datasets
  ) %>%
  relocate(Sample_id:Age_AR, .before = 4)


# readjust site order if needed
merged_data <- merged_data %>%
  mutate(site_order = case_when(
    Site == "Kal0" ~ 1,  
    Site == "Kal1" ~ 2,
    Site == "Kal2" ~ 3,
    Site == "Ran0" ~ 4,
    Site == "Ran1" ~ 5,
    Site == "Ran2" ~ 6,
    Site == "SE17" ~ 7,
    Site == "A5" ~ 8,
    Site == "A13" ~ 9,
    Site == "SE1" ~ 10,
    Site == "Ore1" ~ 11,
    Site == "Ore2" ~ 12,
    Site == "GA1" ~ 13,
    Site == "Gav1" ~ 14,
    Site == "Gav2" ~ 15,
    Site == "C3" ~ 16,
    Site == "C14" ~ 17,
    Site == "SE3" ~ 18,
    Site == "Ljs0" ~ 19,
    Site == "Ljs1" ~ 20,
    Site == "Ljs2" ~ 21,
    Site == "SE4" ~ 22,
    TRUE ~ NA_integer_  # fallback in case of unknown site names
  ))


merged_data <- merged_data %>%
  relocate(site_order, .before = 2)





    #### Adjusting ----

# Here we create a subset of the "merged_data" dataset (without the spectral data) useful for our time series modeling. 
# We start by flooring the age calculated using sedimentation rates (see section 4.1 for more details). 
# Then we add our categories: Distance_to_shore, Basin and Transects (normal and extended).
# Finally, for our spatial variability, we select the top 3 cm, ie. the surface sediment. 

data_by_year <- merged_data[1:68]    # keep relevant columns
data_by_year$Age_AR <- floor(data_by_year$Age_AR)    # convert age to whole years (ie. 1860 instead of 1860.7)

data_by_year <- data_by_year %>%    # add column sorting sites into offshore or coastal 
  mutate(Distance_to_shore = case_when(
    Site %in% c("Kal0", "Kal1", "Kal2",
                "Ran0", "Ran1", "Ran2",
                "Ore1", "Ore2",
                "GA1", "Gav1", "Gav2",
                "Ljs0", "Ljs1", "Ljs2") ~ "Coastal",
    Site %in% c("SE17", "SE1", "SE3", "SE4",
                "A5", "A13", "C3", "C14") ~ "Offshore",
  ))

data_by_year$Distance_to_shore <- as.factor(data_by_year$Distance_to_shore)    # convert to factor


data_by_year <- data_by_year %>%    # add column sorting sites into bay (north) or sea (south)
  mutate(Basin = case_when(
    Site %in% c("Kal0", "Kal1", "Kal2",
                "Ran0", "Ran1", "Ran2",
                "SE17", "SE1",
                "A5", "A13") ~ "Bothnian Bay",
    Site %in% c("Ore1", "Ore2",
                "GA1", "Gav1", "Gav2",
                "Ljs0", "Ljs1", "Ljs2",
                "SE3", "SE4",
                "C3", "C14") ~ "Bothnian Sea",
  ))

data_by_year$Basin <- as.factor(data_by_year$Basin)


data_by_year <- data_by_year %>%    # add column sorting sites into transects
  mutate(Transect = case_when(
    Site %in% c("Kal0", "Kal1", "Kal2") ~ "Kalix", 
    Site %in% c("Ran0", "Ran1", "Ran2") ~ "Ranea",
    Site %in% c("Ore1", "Ore2") ~ "Ore",
    Site %in% c("GA1", "Gav1", "Gav2") ~ "Gavik",
    Site %in% c("Ljs0", "Ljs1", "Ljs2") ~ "Ljusnan",
    Site %in% c("A5", "A13", "SE17", "SE1") ~ "Bothnian Bay",
    Site %in% c("C3", "C14","SE3", "SE4") ~ "Bothnian Sea",
  ))

data_by_year$Transect <- as.factor(data_by_year$Transect)


data_by_year <- data_by_year %>%    # add column sorting sites into extended transects (including offshore points)
  mutate(Extended = case_when(
    Site %in% c("Kal0", "Kal1", "Kal2", "SE17") ~ "Kalix", 
    Site %in% c("Ran0", "Ran1", "Ran2", "A5") ~ "Ranea",
    Site %in% c("Ore1", "Ore2") ~ "Ore",
    Site %in% c("GA1", "Gav1", "Gav2", "C3") ~ "Gavik",
    Site %in% c("Ljs0", "Ljs1", "Ljs2", "SE3") ~ "Ljusnan"
  ))

data_by_year$Extended <- as.factor(data_by_year$Extended)


# here we create a dataset for spatial variability, ie. top 3 cm
# we also adjust the monitoring dataset to inlcude the same groups as in data_by_year

surface_sed <- data_by_year %>%
  filter(Depth <= 3)

monitoring <- monitoring %>%    # add column sorting sites into bay (north) or sea (south)
  mutate(Basin = case_when(
    Site %in% c("Ran2", "A13") ~ "Bothnian Bay",
    Site %in% c("GA1", "B3", "C3") ~ "Bothnian Sea",
  ))

monitoring$Basin <- as.factor(monitoring$Basin)

monitoring <- monitoring %>%    # add column sorting sites into bay (north) or sea (south)
  mutate(Distance_to_shore = case_when(
    Site %in% c("Ran2", "GA1") ~ "Coastal",
    Site %in% c("A13", "B3", "C3") ~ "Offshore",
  ))

monitoring$Distance_to_shore <- as.factor(monitoring$Distance_to_shore)

monitoring <- monitoring %>%
  mutate(site_order = case_when(
    Site == "Ran2" ~ 1,  
    Site == "A13" ~ 2,
    Site == "B3" ~ 3,
    Site == "GA1" ~ 4,
    Site == "C3" ~ 5,
  ))

overlap <- surface_sed %>%
  filter(Site %in% c("Ran2", "GA1", "Ore2", "A13", "C3"))
# overlap allows us to select the same sites as for the monitoring data,
# therefore we can better compare the 2


# for our boxplots (see section 2.2, 2.3, 2.4 and 2.5) we create a colour palette 

surface_sed$Basin_Distance <- interaction(surface_sed$Basin, surface_sed$Distance_to_shore)

colors <- c(
  "Bothnian Bay.Coastal" = "#1b9e77",   
  "Bothnian Bay.Offshore"  = "#70c4b8",   
  "Bothnian Sea.Coastal" = "#7570b3",   
  "Bothnian Sea.Offshore"  = "#b7aef0"    
)

basin_colors <- c(
  "Bothnian Bay" = "#1b9e77",
  "Bothnian Sea" = "#7570b3"
)

distance_colors <- c(
  "Coastal"  = "grey50",
  "Offshore" = "grey80"
)





### 1.4 Quality control ----

# Checking for data consistency (no shift over time)

    #### XRF ----

# loop to go through all elements
for (element in names(XRF_ref)[10:34]) {        #  element is a loop variable, a placeholder that takes on each value from the vector: vector contains the names of all the columns in dataset from column 7 to 17
  p <- ggplot(XRF_ref, aes(x = order, y = .data[[element]])) +
    geom_point() +
    geom_smooth(method = "lm", color = "red") +
    facet_wrap(~ Site) +       # group per site
    labs(
      x = "Approximate run order",
      y = paste0(element, " (dry weight %)"),
      title = paste("Drift check for", element, "across reference materials")
    )
  print(p)
}


# only plotting selected elements
plots <- list()  # empty list to store plots

for (element in c("Mn", "Zn", "Pb")) {
  p <- ggplot(XRF_ref, aes(x = order, y = .data[[element]])) +
    geom_point() +
    geom_smooth(method = "lm", color = "red") +
    facet_wrap(~ Site) +
    labs(
      x = "Approximate run order",
      y = paste0(element, " (dry weight %)"),
      title = paste("Drift check for", element, "across reference materials")
    )
  plots[[element]] <- p   # store each plot in the list
}


plots$Mn / plots$Zn / plots$Pb





    #### VNIR ----

g1 <- ggplot(VNIR_ref, aes(x = order, y = TOC_Uni)) +
  geom_point() +
  geom_smooth(method = "lm", color = "red") +
  facet_wrap(~ Site) +
  labs(
    x = "Approximate run order",
    y = "Inferred WC-TOC",
    title = "Drift check for inferred WC-TOC across reference materials"
  )

# Plot 2
g2 <- ggplot(VNIR_ref, aes(x = order, y = Chla_corr)) +
  geom_point() +
  geom_smooth(method = "lm", color = "red") +
  facet_wrap(~ Site) +
  labs(
    x = "Approximate run order",
    y = "Chl-a",
    title = "Drift check for chl-a across reference materials"
  )

# Combine side-by-side or top-bottom
g1 / g2  # Top-bottom layout
# OR use p1 | p2 for side-by-side





# 2 Spatial variability ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

# For spatial variability, we take the top 3cm from our sediment cores and compare variables in space.
# We look at the difference from coastal to offshore environments, in the different basins 
# and at individual transects/sites, especially with a N to S distribution.
# We also want to compare our surface sediment data to the UMF monitoring data. 

# For total phytoplankton, we use Chla as a biomass indicator, inferred from VNIR spectroscopy
# For diatoms, we use the Si/detritals ratio as an indicator, from XRF analysis
# For comparison with monitoring data, we select the same stations as for the monitoring data 
# (ie. Ran2, A13, GA1, Ore2, C3) and compare chla from sediments with total phytoplankton biomass from UMF. 



### 2.1 Data exploration ----

# We start by exploring the data using a PCA to see whether we can detect any patterns

surfaceList <- surface_sed %>%      # creating a list of diagenetic elements
  select(TOC_Uni, Chla_corr, Si_detrital, P, Mn_Fe, S) %>%  # keep relevant columns
  names()

# for the PCA we use scaled and centered data transformation because the variables don't have the same units

# plotting PC1 vs PC2
p1 <- surface_sed %>%
  arrange(site_order) %>%  # First, order rows by your custom numeric order
  mutate(Site = factor(Site, levels = unique(Site))) %>%  # Then turn Site into an ordered factor
  ordr::ordinate(
    cols = any_of(surfaceList),
    model = ~ prcomp(., scale. = TRUE), 
    #model = ~ princomp(x = ., cor = TRUE, scores = TRUE),  
    augment = any_of(c("Depth", "Site", "Basin"))
  ) %>%
  ordr::ggbiplot(., mapping = aes(x = 1, y = 2), sec.axes = "cols", scale.factor = 6 #, groups = CONISS, ellipse = TRUE, ellipse.prob = 0.95, circle = TRUE
  ) +
  ordr::geom_rows_point(aes(colour = Site, shape = Basin), alpha = 0.7, size = 2) +
  ordr::geom_cols_vector() +
  ordr::geom_cols_text_radiate(aes(label = name), size = 3) +
  scale_color_viridis_d(option = "C") +
  #scale_color_hue() +
  labs(title = "PCA of Surface Sediments", color = "Site", shape = "Basin") +
  xlim(-5,5) +
  ylim(-5,5) +
  theme_minimal() +
  theme(legend.position = "bottom")

# plotting PC2 vs PC3
p2 <- surface_sed %>%
  arrange(site_order) %>%  # Order rows by custom numeric order
  mutate(Site = factor(Site, levels = unique(Site))) %>%  # Ordered factor for Site
  ordr::ordinate(
    cols = any_of(surfaceList),
    #model = ~ prcomp(clr(.)),
    model = ~ prcomp(., scale. = TRUE), 
    augment = any_of(c("Depth", "Site", "Basin"))
  ) %>%
  ordr::ggbiplot(mapping = aes(x = 2, y = 3), sec.axes = "cols", scale.factor = 3) +  # PC2 vs PC3
  ordr::geom_rows_point(aes(colour = Site, shape = Basin), alpha = 0.7, size = 2) +
  ordr::geom_cols_vector() +
  ordr::geom_cols_text_radiate(aes(label = name), size = 3) +
  scale_color_viridis_d(option = "C") +
  #labs(title = "PCA of Surface Sediments (PC2 vs PC3)",
  #color = "Site", shape = "Basin") +
  xlim(-3,3) +
  ylim(-4,3) +
  theme_minimal() +
  theme(legend.position = "right")


p1 + p2 

# we notice that PCA of surface sediments shows clustering of TOC and Mn/Fe towards the Bothnian Bay, 
# whereas Chl-a looks to be pointing towards the Bothnian Sea




### 2.2 Monitoring vs phytoplankton ----

# here we plot boxplots of monitoring phytoplankton biomass vs our chl-a (which we use as a phytoplankton
# abundance indicator) from surface sediments, comparing both basins and using the sites that overlap 
# (ie. the sites from the sediment data that are the same as the sites from the monitoring data)

q1 <- ggplot(monitoring, aes(x = Basin, y = Phytoplankton, fill = Basin)) +
  geom_boxplot(alpha = 0.7) +
  scale_fill_manual(values = basin_colors) +
  theme_bw() +
  labs(
    title = "Monitoring phytoplankton biomass",
    x = "Basin",
    y = "Weighted yearly means of 
  total phytoplankton (mg/m3)"
  ) +
  theme(legend.position = "none")

q2 <- ggplot(overlap, aes(x = Basin, y = Chla_corr, fill = Basin)) +
  geom_boxplot(alpha = 0.7) +
  scale_fill_manual(values = basin_colors) +
  theme_bw() +
  labs(
    title = "Chl-a inferred from surface sediments",
    x = "Basin",
    y = "Chl-a concentration (mg/g dry mass)"
  ) +
  theme(legend.position = "none")

q1 + q2


# to test the difference between basins for each variable we use a Wilcoxon rank sum test

wilcox.test(Phytoplankton ~ Basin, data = monitoring)

monitoring %>%
  group_by(Basin) %>%
  summarise(
    median_phyto = median(Phytoplankton, na.rm = TRUE),
    IQR_phyto = IQR(Phytoplankton, na.rm = TRUE),
    .groups = "drop"
  )

wilcox.test(Chla_corr ~ Basin, data = overlap)

overlap %>%
  group_by(Basin) %>%
  summarise(
    median_chla = median(Chla_corr, na.rm = TRUE),
    IQR_chla = IQR(Chla_corr, na.rm = TRUE),
    .groups = "drop"
  )





### 2.3 Phytoplankton ----

# Here we look at the data from all of our sites, not only the ones that coincide with the monitoring data

    #### Total Chla ----

# for total phytoplankton abundance, we use chl-a as an indicator, inferred from VNIR spectroscopy

ggplot(surface_sed, aes(x = "", y = Chla_corr, fill = Basin_Distance)) +  # x is a constant: one box per panel
  geom_boxplot(alpha = 0.7) +
  facet_grid(Basin ~ Distance_to_shore) +  # 2x2 panels
  scale_fill_manual(values = colors) +
  theme_bw() +
  theme(axis.text.x = element_blank(),  # no x-axis labels
        axis.ticks.x = element_blank()) +
  labs(title = "Chl-a inferred from surface sediments",
       x = NULL,
       y = "Chl-a concentration (mg/g dry mass)")


# to test the difference between groups we use a kruskal wallis test followed by a pairwise test

kruskal_chla <- kruskal_test(Chla_corr ~ Basin_Distance, data = surface_sed)
print(kruskal_chla)

pairwise_chla <- surface_sed %>%
  pairwise_wilcox_test(Chla_corr ~ Basin_Distance, p.adjust.method = "holm")
print(pairwise_chla)





    #### Diatoms ----

# for diatoms we use the Si/detritals ratio as an indicator, calculated from XRF scanning

ggplot(surface_sed, aes(x = "", y = Si_detrital, fill = Basin_Distance)) +  # x is a constant: one box per panel
  geom_boxplot(alpha = 0.7) +
  facet_grid(Basin ~ Distance_to_shore) +  # 2x2 panels
  scale_fill_manual(values = colors) +
  theme_bw() +
  theme(axis.text.x = element_blank(),  # no x-axis labels
        axis.ticks.x = element_blank()) +
  labs(title = "Diatoms inferred from surface sediments",
       x = NULL,
       y = "Si/detritals")


# to test the difference between groups we use a kruskal wallis test followed by a pairwise test

kruskal_diat <- kruskal_test(Si_detrital ~ Basin_Distance, data = surface_sed)
print(kruskal_diat)

pairwise_diat <- surface_sed %>%
  pairwise_wilcox_test(Si_detrital ~ Basin_Distance, p.adjust.method = "holm")
print(pairwise_diat)





### 2.4 Water colour ----

# For water colour, we use water column TOC as an indicator, inferred from VNIR spectroscopy
# there being not coastal-offshore gradient, we plot only the difference between basins

ggplot(surface_sed, aes(x = Basin, y = TOC_Uni, fill = Basin)) +
  geom_boxplot(alpha = 0.7) +
  theme_bw() +
  labs(
    title = "WC-TOC inferred from surface sediments",
    x = "Basin",
    y = "WC-TOC (mg/L)"
  ) +
  scale_fill_manual(values = basin_colors) +
  theme(legend.position = "none")


# to test the difference between the 2 basins we use a t-test, since the preliminary conditions are met

# normality
with(surface_sed, shapiro.test(TOC_Uni[Basin == "Bothnian Bay"]))
with(surface_sed, shapiro.test(TOC_Uni[Basin == "Bothnian Sea"]))

ggplot(surface_sed, aes(sample = TOC_Uni)) +
  stat_qq() + stat_qq_line() +
  facet_wrap(~ Basin)

# homogeneity of variances
leveneTest(TOC_Uni ~ Basin, data = surface_sed)

# t-test
t.test(TOC_Uni ~ Basin, data = surface_sed, var.equal = TRUE)





### 2.5 Redox conditions ----

# for oxygen levels we use the Mn/Fe ratio as an indicator, calculated from XRF scanning
# low ratios indicate low oxygen levels and vice versa

ggplot(surface_sed, aes(x = "", y = Mn_Fe, fill = Basin_Distance)) +  # map fill to combined factor
  geom_boxplot(alpha = 0.7) +
  facet_grid(Basin ~ Distance_to_shore) +  # 2x2 panels
  scale_fill_manual(values = colors) +
  theme_bw() +
  theme(axis.text.x = element_blank(),  # no x-axis labels
        axis.ticks.x = element_blank()) +
  labs(title = "Redox conditions inferred from surface sediments",
       x = NULL,
       y = "Mn/Fe")


# to test the difference between groups we use a kruskal wallis test followed by a pairwise test

kruskal_redox <- kruskal_test(Mn_Fe ~ Basin_Distance, data = surface_sed)
print(kruskal_redox)

pairwise_redox <- surface_sed %>%
  pairwise_wilcox_test(Mn_Fe ~ Basin_Distance, p.adjust.method = "holm")
print(pairwise_redox)





# 3 Temporal variability ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

# For temporal variability, we take all the data from our sediment cores and compare variables in time
# We look at the difference from coastal to offshore environments and in the different basins.
# For plotting, the default "plot()" function is usable to plot gams.
# However, to make aesthetic plots (customised, more pretty, etc.), we extract the data from the 
# gam model and use ggplot() - see below. 


### 3.1 Phytoplankton ----

# For total phytoplankton, we use Chla as a biomass indicator, inferred from VNIR spectroscopy
# For diatoms, we use the Si/detritals ratio as an indicator, from XRF analysis


    #### Total Chla ----

# We start with the global trend of Chla, ie. using the data from all cores in the Gulf of Bothnia combined.
# Then we look at groupings (offshore vs coastal, bay vs sea). 


############ global trend

global_chla <- data_by_year[, c("Age_AR", "Chla_corr")]

# average data points per year to reduce noise
global_chla <- global_chla %>%     
  group_by(Age_AR) %>%
  summarise(mean_chla = mean(Chla_corr, na.rm = TRUE))

# run model
gam_global_chla <- gam(mean_chla ~ s(Age_AR, k=9, bs='tp'),
                       data =  subset(global_chla, Age_AR >= 1850), 
                       method = "REML")

summary(gam_global_chla)

# add gam data to dataset
chla_gl <- subset(global_chla, Age_AR >= 1850)
chla_gl$fit <- predict(gam_global_chla, type = "response")
chla_gl$se <- predict(gam_global_chla, type = "response", se.fit = TRUE)$se.fit

# plot using 
ggplot(chla_gl, aes(x = Age_AR, y = fit)) +
  # Confidence ribbon
  geom_ribbon(aes(ymin = fit - 2 * se, ymax = fit + 2 * se), 
              fill = "paleturquoise3", alpha = 0.4) +
  # Fitted line
  geom_line(color = "paleturquoise4", size = 1.2) +
  # Rug plot to show data density
  geom_rug(aes(x = Age_AR), sides = "b", alpha = 0.3, color = "gray40", length = unit(0.02, "npc")) +
  # Optional: add actual data points
  geom_point(aes(y = mean_chla), color = "black", size = 1, alpha = 0.5) +
  # Vertical line
  #geom_vline(xintercept = 1920, linetype = "dashed", color = "red", size = 1) +
  # Labels
  labs(
    title = "Chl-a trend over time in the Gulf of Bothnia",
    x = "Age (years)",
    y = "Mean Chl-a (mg/g dry mass)"
  ) +
  # Theme improvements
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 15),
    axis.title = element_text(size = 13, face = "bold"),
    panel.grid.minor = element_blank()
  ) 


############ coastal vs offshore

gam_dist_chla <- gam(Chla_corr ~ s(Age_AR, by = Distance_to_shore, k = 9, bs ="fs"), 
                     data = subset(data_by_year, Age_AR >= 1850),
                     method = "REML")

summary(gam_dist_chla)

# to plot, simply use:
plot(gam_dist_chla, pages = 1, residuals = TRUE)
# or use code for "pretty" plot as above 


############ bay vs sea

gam_basin_chla <- gam(Chla_corr ~ s(Age_AR, by = Basin, k = 9, bs ="fs"), 
                      data = subset(data_by_year, Age_AR >= 1850),
                      method = "REML")

summary(gam_basin_chla)

# to plot, simply use:
plot(gam_basin_chla, pages = 1, residuals = TRUE)
# or use code for "pretty" plot as above 





# # # # # # # # # # # 
# Adapt code depending on what you want to look at.
# # # # # # # # # # # 





    #### Diatoms ----

# Same procedure as with Chla. 


############ global trend

global_diat <- data_by_year[, c("Age_AR", "Si_detrital")]

# average data points per year to reduce noise
global_diat <- global_diat %>%
  group_by(Age_AR) %>%
  summarise(mean_diat = mean(Si_detrital, na.rm = TRUE))

# run model
gam_global_diat <- gam(mean_diat ~ s(Age_AR, k=9, bs='tp'),
                       data =  subset(global_diat, Age_AR >= 1850), 
                       method = "REML")

summary(gam_global_diat)

plot(gam_global_diat, residuals = TRUE, main = paste("Si/detritals time series for all sites"))


############ coastal vs offshore

gam_dist_diat <- gam(Si_detrital ~ s(Age_AR, by = Distance_to_shore, k = 9, bs ="fs"), 
                     data = subset(data_by_year, Age_AR >= 1850),
                     method = "REML")

summary(gam_dist_diat)

plot(gam_dist_diat, residuals = TRUE, main = paste("Si/detritals time series for all sites"))



############ bay vs sea

gam_basin_diat <- gam(Si_detrital ~ s(Age_AR, by = Basin, k = 9, bs ="fs"), 
                      data = subset(data_by_year, Age_AR >= 1850),
                      method = "REML")

summary(gam_basin_diat)

diat_basin <- subset(data_by_year, Age_AR >= 1850)

# Add gam data to dataset
pred1 <- predict(gam_basin_diat, newdata = diat_basin, type = "response", se.fit = TRUE)
diat_basin$fit <- pred1$fit
diat_basin$se <- pred1$se.fit

# Plot with color per basin
ggplot(diat_basin, aes(x = Age_AR, y = fit, color = Basin, fill = Basin)) +
  # Confidence ribbon
  geom_ribbon(aes(ymin = fit - 2 * se, ymax = fit + 2 * se), alpha = 0.3, color = NA) +
  # Fitted line
  geom_line(size = 1.2) +
  # Rug plot
  geom_rug(aes(x = Age_AR), sides = "b", alpha = 0.3, color = "gray40", length = unit(0.02, "npc")) +
  # Raw points (optional)
  # geom_point(aes(y = Si_detrital), color = "black", size = 1, alpha = 0.5) +
  # Vertical line
  #geom_vline(xintercept = 1920, linetype = "dashed", color = "red", size = 1) +
  #geom_vline(xintercept = 1975, linetype = "dashed", color = "blue", size = 1) +
  # Labels
  labs(
    title = "Diatom trend over time in the Gulf of Bothnia",
    x = "Age (years)",
    y = "Si/detritals",
    color = "Basin",
    fill = "Basin"
  ) +
  # Color palette
  scale_color_manual(values = c("Bothnian Bay" = "#1b9e77", "Bothnian Sea" = "#7570b3")) +
  scale_fill_manual(values = c("Bothnian Bay" = "#1b9e77", "Bothnian Sea" = "#7570b3")) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 15),
    axis.title = element_text(size = 13, face = "bold"),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )


# # # # # # # # # # # 
# Here, the trend is different between bay and sea, which explains why the overall model 
# is not very significant/ deviance explained is not very high. 
# We can test whether the coastal vs offshore response is different in one of the basins. 
# # # # # # # # # # # 


############ bay: coastal vs offshore

gam_bay_diat <- gam(Si_detrital ~ s(Age_AR, by = Distance_to_shore, k = 9, bs ="fs"), 
                    data = subset(data_by_year, Basin == "Bothnian Bay" & Age_AR >= 1850),
                    method = "REML")

summary(gam_bay_diat)

# Ensure Distance_to_shore is a factor
diat_bay <- subset(data_by_year, Basin == "Bothnian Bay" & Age_AR >= 1850)
diat_bay$Distance_to_shore <- as.factor(diat_bay$Distance_to_shore)

# Add gam data to dataset
pred <- predict(gam_bay_diat, newdata = diat_bay, type = "response", se.fit = TRUE)
diat_bay$fit <- pred$fit
diat_bay$se <- pred$se.fit

# Plot with color per Distance_to_shore
ggplot(diat_bay, aes(x = Age_AR, y = fit, color = Distance_to_shore, fill = Distance_to_shore)) +
  # Confidence ribbon
  geom_ribbon(aes(ymin = fit - 2 * se, ymax = fit + 2 * se), alpha = 0.3, color = NA) +
  # Fitted line
  geom_line(size = 1.2) +
  # Rug plot
  geom_rug(aes(x = Age_AR), sides = "b", alpha = 0.3, color = "gray40", length = unit(0.02, "npc")) +
  # Raw points (optional)
  # geom_point(aes(y = Si_detrital), color = "black", size = 1, alpha = 0.5) +
  # Labels
  labs(
    title = "Si/detritals trend over time in the Bothnian Bay",
    x = "Age (years)",
    y = "Si/detritals",
    color = "Distance to shore",
    fill = "Distance to shore"
  ) +
  # Color palette
  scale_color_manual(values = c("Coastal" = "#1b9e77", "Offshore" = "#7570b3")) +
  scale_fill_manual(values = c("Coastal" = "#1b9e77", "Offshore" = "#7570b3")) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 15),
    axis.title = element_text(size = 13, face = "bold"),
    panel.grid.minor = element_blank()
  )





### 3.2 Browning ----

# For browning, we use water column TOC as an indicator, inferred from VNIR spectroscopy

# # # #   /!\   # # # # 
#
# Here, we exclude the top layer in sediments subject to the most intensive post-depositional 
# alteration, diagenetic remobilization and microbial activity. 
# This top layer is around 3cm deep and equivalent to the last ~20 years.
# Both TOC and diagenetic elements, but not chla, appeared most sensitive to these processes, 
# as visible from the ‘hockey-stick’ shift in trends downcore. 
# More in section 4.4
# 
# # # #   /!\   # # # # 


############ global trend

global_toc <- data_by_year[, c("Age_AR", "TOC_Uni")]

global_toc <- global_toc %>%
  group_by(Age_AR) %>%
  summarise(mean_toc = mean(TOC_Uni, na.rm = TRUE))


gam_global_toc <- gam(mean_toc ~ s(Age_AR, k=9, bs='tp'),
                      data =  subset(global_toc, Age_AR >= 1850), # cut off last 20y = diagenetic layer
                      method = "REML")

plot(gam_global_toc, residuals = TRUE, main = paste("TOC time series for all sites"))


summary(gam_global_toc)

toc_gl <- subset(global_toc, Age_AR >= 1850)
toc_gl$fit <- predict(gam_global_toc, type = "response")
toc_gl$se <- predict(gam_global_toc, type = "response", se.fit = TRUE)$se.fit

# Plot

ggplot(toc_gl, aes(x = Age_AR, y = fit)) +
  # Confidence ribbon
  geom_ribbon(aes(ymin = fit - 2 * se, ymax = fit + 2 * se), 
              fill = "paleturquoise3", alpha = 0.4) +
  # Fitted line
  geom_line(color = "paleturquoise4", size = 1.2) +
  # Rug plot to show data density
  geom_rug(aes(x = Age_AR), sides = "b", alpha = 0.3, color = "gray40", length = unit(0.02, "npc")) +
  # Optional: add actual data points
  geom_point(aes(y = mean_toc), color = "black", size = 1, alpha = 0.5) +
  # 🔴 Cross-out area
  annotate("rect",
           xmin = max(toc_gl$Age_AR) - 23,
           xmax = max(toc_gl$Age_AR),
           ymin = -Inf, ymax = Inf,
           fill = "red", alpha = 0.2) +
  # Labels
  labs(
    title = "     WC-TOC trend over time in the Gulf of Bothnia",
    x = "Age (years)",
    y = "Mean inferred WC-TOC (mg/L)"
  ) +
  # Theme improvements
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 15),
    axis.title = element_text(size = 13, face = "bold"),
    panel.grid.minor = element_blank()
  ) + 
  ylim(2,5)


############ coastal vs offshore

gam_dist_toc <- gam(TOC_Uni ~ s(Age_AR, by = Distance_to_shore, k = 9, bs ="fs"), 
                    data = subset(data_by_year, Age_AR >= 1850 & Age_AR <= 2000),
                    method = "REML")

plot(gam_dist_toc, pages = 1, residuals = TRUE)

summary(gam_dist_toc)


############ bay vs sea

gam_basin_toc <- gam(TOC_Uni ~ s(Age_AR, by = Basin, k = 9, bs ="fs"), 
                     data = subset(data_by_year, Age_AR >= 1850 & Age_AR <= 2000),
                     method = "REML")

plot(gam_basin_toc, pages = 1, residuals = TRUE)

summary(gam_basin_toc)


# # # # # # # # # # # 
# Here, the only time series with a significant trend is the offshore one.
# Let's plot this one separately.
# # # # # # # # # # # 


############ offshore only

gam_offshore_toc <- gam(TOC_Uni ~ s(Age_AR, k = 9, bs ="fs"), 
                        data = subset(data_by_year, Distance_to_shore == "Offshore" & Age_AR >= 1850 & Age_AR <= 2000),
                        method = "REML")

summary(gam_offshore_toc)

# Add gam data to dataset
toc_off <- subset(data_by_year, Distance_to_shore == "Offshore" & Age_AR >= 1850 & Age_AR <= 2000)
toc_off$fit <- predict(gam_offshore_toc, type = "response")
toc_off$se <- predict(gam_offshore_toc, type = "response", se.fit = TRUE)$se.fit

# Plot

ggplot(toc_off, aes(x = Age_AR, y = fit)) +
  # Confidence ribbon
  geom_ribbon(aes(ymin = fit - 2 * se, ymax = fit + 2 * se), 
              fill = "grey80", alpha = 0.4) +
  # Fitted line
  geom_line(color = "grey40", size = 1.2) +
  # Rug plot to show data density
  geom_rug(aes(x = Age_AR), sides = "b", alpha = 0.3, color = "gray40", length = unit(0.02, "npc")) +
  # Optional: add actual data points
  geom_point(aes(y = TOC_Uni), color = "black", size = 1, alpha = 0.5) +
  # Vertical line
  #geom_vline(xintercept = 1940, linetype = "dashed", color = "red", size = 1) +
  # Labels
  labs(
    title = "     WC-TOC trend over time in outer-sea Gulf of Bothnia",
    x = "Age (years)",
    y = "Inferred WC-TOC (mg/L)"
  ) +
  # Theme improvements
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 15),
    axis.title = element_text(size = 13, face = "bold"),
    panel.grid.minor = element_blank()
  ) + 
  ylim(2,5)





### 3.3 Redox conditions ----

############ global trend

global_redox <- data_by_year[, c("Age_AR", "Mn_Fe")]

global_redox <- global_redox %>%
  group_by(Age_AR) %>%
  summarise(mean_redox = mean(Mn_Fe, na.rm = TRUE))


gam_global_redox <- gam(mean_redox ~ s(Age_AR, k=9, bs='tp'),
                        data =  subset(global_redox, Age_AR >= 1850 & Age_AR <= 2000), 
                        method = "REML")

summary(gam_global_redox)

# Add gam data to dataset
redox_gl <- subset(global_redox, Age_AR >= 1850 & Age_AR <= 2000)
redox_gl$fit <- predict(gam_global_redox, type = "response")
redox_gl$se <- predict(gam_global_redox, type = "response", se.fit = TRUE)$se.fit

# Plot

ggplot(redox_gl, aes(x = Age_AR, y = fit)) +
  # Confidence ribbon
  geom_ribbon(aes(ymin = fit - 2 * se, ymax = fit + 2 * se), 
              fill = "paleturquoise3", alpha = 0.4) +
  # Fitted line
  geom_line(color = "paleturquoise4", size = 1.2) +
  # Rug plot to show data density
  geom_rug(aes(x = Age_AR), sides = "b", alpha = 0.3, color = "gray40", length = unit(0.02, "npc")) +
  # Optional: add actual data points
  geom_point(aes(y = mean_redox), color = "black", size = 1, alpha = 0.5) +
  # Labels
  labs(
    title = "Global Mn/Fe trend over time in the Gulf of Bothnia",
    x = "Age (years)",
    y = "Mn/Fe"
  ) +
  # Theme improvements
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 15),
    axis.title = element_text(size = 13, face = "bold"),
    panel.grid.minor = element_blank()
  )


############ coastal vs offshore

gam_dist_redox <- gam(Mn_Fe ~ s(Age_AR, by = Distance_to_shore, k = 9, bs ="fs"), 
                      data = subset(data_by_year, Age_AR >= 1850 & Age_AR <= 2000),
                      method = "REML")

plot(gam_dist_redox, pages = 1, residuals = TRUE)

summary(gam_dist_redox)


############ bay vs sea

gam_basin_redox <- gam(Mn_Fe ~ s(Age_AR, by = Basin, k = 9, bs ="fs"), 
                       data = subset(data_by_year, Age_AR >= 1850 & Age_AR <= 2000),
                       method = "REML")

plot(gam_basin_redox, pages = 1, residuals = TRUE)

summary(gam_basin_redox)


# # # # # # # # # # # 
# Here, the global model is significant but not by much, with a deviance explained of 6.06%.
# The coastal vs offshore model is not significant. The Bothnian Bay model is similar to the global one.  
# If we plot the Bothnian Bay by distance_to_shore, the offshore one seems significant.
# Let's plot this one separately.
# # # # # # # # # # # 


############ offshore bay

gam_bay_redox <- gam(Mn_Fe ~ s(Age_AR, k = 9, bs ="fs"), 
                     data = subset(data_by_year, Basin == "Bothnian Bay" & Distance_to_shore == "Offshore" 
                                   & Age_AR >= 1850 & Age_AR <= 2000),
                     method = "REML")

summary(gam_bay_redox)

# Add gam data to dataset
redox_bay <- subset(data_by_year, Basin == "Bothnian Bay" & Distance_to_shore == "Offshore" 
                    & Age_AR >= 1850 & Age_AR <= 2000)
redox_bay$fit <- predict(gam_bay_redox, type = "response")
redox_bay$se <- predict(gam_bay_redox, type = "response", se.fit = TRUE)$se.fit

# Plot

ggplot(redox_bay, aes(x = Age_AR, y = fit)) +
  # Confidence ribbon
  geom_ribbon(aes(ymin = fit - 2 * se, ymax = fit + 2 * se), 
              fill = "paleturquoise3", alpha = 0.4) +
  # Fitted line
  geom_line(color = "paleturquoise4", size = 1.2) +
  # Rug plot to show data density
  geom_rug(aes(x = Age_AR), sides = "b", alpha = 0.3, color = "gray40", length = unit(0.02, "npc")) +
  # Optional: add actual data points
  geom_point(aes(y = Mn_Fe), color = "black", size = 1, alpha = 0.5) +
  # Labels
  labs(
    title = "Mn/Fe trend over time in outer-sea Bothnian Bay",
    x = "Age (years)",
    y = "Mn/Fe"
  ) +
  # Theme improvements
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 15),
    axis.title = element_text(size = 13, face = "bold"),
    panel.grid.minor = element_blank()
  )



# 4 Modeling ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

# Here we want to see whether TOC drives phytoplankton in the Gulf of Bothnia.
# Similarly as before, we use GAM to model Chla (for total phytoplankton) against TOC, 
# then Si/detritals (for diatoms) against TOC. 
# We also look at groupings such as coastal vs offshore, bay vs sea.
# Don't forget to cut off that diagenetic layer for TOC to avoid misinterpretation. 

### 4.1 TOC vs total chla ----

############ global trend

model_global <- gam(Chla_corr ~ s(TOC_Uni), data = subset(merged_data, TOC_Uni <= 4 & Age_AR >= 1850 & Age_AR <= 2000))


plot(model_global, pages = 1, residuals = TRUE)

summary(model_global)


############ coastal vs offshore

model_dist <- gam(Chla_corr ~ s(TOC_Uni, by = Distance_to_shore, k = 9, bs ="fs"), 
                  data = subset(data_by_year, TOC_Uni <= 4.5 & Age_AR >= 1850 & Age_AR <= 2000))

plot(model_dist, pages = 1, residuals = TRUE)

summary(model_dist)


############ bay vs sea

model_basin <- gam(Chla_corr ~ s(TOC_Uni, by = Basin, k = 9, bs ="fs"), 
                   data = subset(data_by_year, TOC_Uni <= 4.5 & Age_AR >= 1850 & Age_AR <= 2000))

plot(model_basin, pages = 1, residuals = TRUE)

summary(model_basin)


# # # # # # # # # # # 
# Here, the only model with a significant trend is the bothnian bay one.
# Let's plot this one separately.
# # # # # # # # # # # 


############ bay only

model_bb <- gam(Chla_corr ~ s(TOC_Uni, k = 9), 
                data = subset(data_by_year, Basin == "Bothnian Bay" & TOC_Uni <= 4 & TOC_Uni >= 2
                              & Age_AR >= 1850 & Age_AR <= 2000))

# Add predictions and SE to your dataset
toc_chla <- subset(data_by_year, Basin == "Bothnian Bay" & TOC_Uni <= 4 & TOC_Uni >= 2 
                   & Age_AR >= 1850 & Age_AR <= 2000)
toc_chla$fit <- predict(model_bb, type = "response")
toc_chla$se <- predict(model_bb, type = "response", se.fit = TRUE)$se.fit

# Plot

ggplot(toc_chla, aes(x = TOC_Uni, y = fit)) +
  # Confidence ribbon
  geom_ribbon(aes(ymin = fit - 2 * se, ymax = fit + 2 * se), 
              fill = "#1b9e77", alpha = 0.2) +
  # Fitted line
  geom_line(color = "#1b9e77", size = 1.2) +
  # Rug plot to show data density
  geom_rug(aes(x = TOC_Uni), sides = "b", alpha = 0.3, color = "gray40", length = unit(0.02, "npc")) +
  # Optional: add actual data points
  geom_point(aes(y = Chla_corr), color = "black", size = 1, alpha = 0.5) +
  # Labels
  labs(
    title = "Model of WC-TOC and chl-a in the Bothnian Bay",
    x = "Inferred WC-TOC (mg/L)",
    y = "Chl-a (mg/g dry mass)"
  ) +
  # Theme improvements
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 15),
    axis.title = element_text(size = 13, face = "bold"),
    panel.grid.minor = element_blank()
  ) 





### 4.1 TOC vs diatoms ----

############ global trend

model_global1 <- gam(Si_detrital ~ s(TOC_Uni), data = subset(data_by_year, TOC_Uni <= 4 & TOC_Uni >= 2.5
                                                             & Age_AR >= 1850 & Age_AR <= 2000))

# Add gam data to dataset
toc_diat <- subset(data_by_year, TOC_Uni <= 4 & TOC_Uni >= 2.5 
                   & Age_AR >= 1850 & Age_AR <= 2000)
toc_diat$fit <- predict(model_global1, type = "response")
toc_diat$se <- predict(model_global1, type = "response", se.fit = TRUE)$se.fit

# Plot

ggplot(toc_diat, aes(x = TOC_Uni, y = fit)) +
  # Confidence ribbon
  geom_ribbon(aes(ymin = fit - 2 * se, ymax = fit + 2 * se), 
              fill = "paleturquoise3", alpha = 0.4) +
  # Fitted line
  geom_line(color = "paleturquoise4", size = 1.2) +
  # Rug plot to show data density
  geom_rug(aes(x = TOC_Uni), sides = "b", alpha = 0.3, color = "gray40", length = unit(0.02, "npc")) +
  # Optional: add actual data points
  geom_point(aes(y = Si_detrital), color = "black", size = 1, alpha = 0.5) +
  # Vertical line
  #geom_vline(xintercept = 3.6, linetype = "dashed", color = "red", size = 1) +
  # Horizontal line
  #geom_hline(yintercept = 5.4, linetype = "dashed", color = "blue", size = 1) +
  # Labels
  labs(
    title = "      Model of WC-TOC and Si/detritals in the Gulf of Bothnia",
    x = "Inferred WC-TOC (mg/L)",
    y = "Si/detritals"
  ) +
  # Theme improvements
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 15),
    axis.title = element_text(size = 13, face = "bold"),
    panel.grid.minor = element_blank()
  ) 

############ coastal vs offshore

model_dist1 <- gam(Si_detrital ~ s(TOC_Uni, by = Distance_to_shore, k = 9), 
                   data = subset(data_by_year, TOC_Uni <= 4.5
                                 & Age_AR >= 1850 & Age_AR <= 2000))

plot(model_dist1, pages = 1, residuals = TRUE)

summary(model_dist1)


############ bay vs sea

model_basin1 <- gam(Si_detrital ~ s(TOC_Uni, by = Basin, k = 9, bs ="fs"), 
                    data = subset(data_by_year, TOC_Uni <= 4.5 
                                  & Age_AR >= 1850 & Age_AR <= 2000))

plot(model_basin1, pages = 1, residuals = TRUE)

summary(model_basin1)





# 5 Misc ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

### 5.1 Age-depth profiles ----

# While waiting for 210Pb dating, we have calculated age using sedimentation rates (EMODnet).
# These profiles show the length of the cores and the age associated with their depth.

ggplot(merged_data, aes(x = Age_AR, y = Depth, color = reorder(Site, site_order))) +
  geom_line(aes(group = Site)) +
  geom_point(size = 2) +
  scale_y_reverse() +
  scale_color_viridis_d(option = "C") + # Options: "A", "B", "C", "D", "E"
  labs(
    title = "Age-Depth profiles by site",
    x = "Age (years)",
    y = "Depth (cm)",
    color = "Site"
  ) +
  theme_minimal()





### 5.2 Lead time series ----

# The lead (Pb) time series reinforces our age calculations, showing the peak in the 1970s,
# which corresponds to the lead pollution in Sweden, followed by measures to reduce emissions. 

global_pb <- data_by_year[, c("Age_AR", "Pb")]

global_pb <- global_pb %>%
  group_by(Age_AR) %>%
  summarise(mean_pb = mean(Pb, na.rm = TRUE)) # averaging reduces noise

gam_global_pb <- gam(mean_pb ~ s(Age_AR, k=9, bs='tp'),
                     data =  subset(global_pb, Age_AR >= 1850), 
                     method = "REML")

summary(gam_global_pb)


# Add gam data to dataset
pb_gl <- subset(global_pb, Age_AR >= 1850)
pb_gl$fit <- predict(gam_global_pb, type = "response")
pb_gl$se <- predict(gam_global_pb, type = "response", se.fit = TRUE)$se.fit


# Plot
ggplot(pb_gl, aes(x = Age_AR, y = fit)) +
  # Confidence ribbon
  geom_ribbon(aes(ymin = fit - 2 * se, ymax = fit + 2 * se), 
              fill = "paleturquoise3", alpha = 0.4) +
  # Fitted line
  geom_line(color = "paleturquoise4", size = 1.2) +
  # Rug plot to show data density
  geom_rug(aes(x = Age_AR), sides = "b", alpha = 0.3, color = "gray40", length = unit(0.02, "npc")) +
  # Optional: add actual data points
  geom_point(aes(y = mean_pb), color = "black", size = 1, alpha = 0.5) +
  # Labels
  labs(
    title = "Global lead (Pb) trend over time in the Gulf of Bothnia",
    x = "Age (years)",
    y = "Mean Pb proportion in sediments (%)"
  ) +
  # Theme improvements
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 15),
    axis.title = element_text(size = 13, face = "bold"),
    panel.grid.minor = element_blank()
  ) 




### 5.3 Testing universal model (TOC) ----

# The universal model for inferring water column TOC from VNIR spectroscopy was developed 
# by Meyer-Jacob et al., 2017 from northern lake sediments. 
# We apply it here for the first time in a (semi-)marine environment.
# To show that it works, we perform a PCA and plot the scores of both spectra. 
# The overlap of our data with the model shows that the model works also for our data. 

spectraList <- VNIR_avg %>%      # creating a list of all wavelengths
  select(where(is.numeric)) %>%
  select(-Depth, -site_order, -TOC_Uni, -TOC_Swe, -VRS_chla, -Chla_corr) %>%  # remove irrelevant columns
  names()


UniModel <- UniModel %>%
  mutate(Site = "UniModel") # adding site column 


VNIR_Unimodel <- bind_rows( # merging both datasets
  VNIR_avg,
  UniModel %>%
    rename(TOC_Uni = TOC)  # Rename to match VNIR_avg
)


VNIR_Unimodel <- VNIR_Unimodel %>%  # add "type" column which indicates where the data comes from 
  mutate(
    Type = case_when(
      is.na(Model)        ~ "Sweden",  # NA becomes Sweden
      Model == "Sweden"   ~ "Sweden",
      Model == "Canada"   ~ "Canada",
      Model == "Finland"  ~ "Finland",
      Model == "Other"    ~ "Other"    
    )
  )


# manually setting colour scale so that the sites (ordered from N to S) follow a colour scale, and the Universal Model is in black
site_colours <- {
  site_levels <- VNIR_Unimodel %>% # make sure the colour scale follows the correct site order
    arrange(site_order) %>%
    pull(Site) %>%
    unique()
  
  setNames(
    ifelse(site_levels == "UniModel", "black", viridis(length(site_levels) - 1, option = "C")),
    site_levels
  )
}


# Create PCA plot with only colour coding
VNIR_Unimodel %>%
  arrange(site_order) %>%  # order rows by site order
  mutate(Site = factor(Site, levels = unique(Site))) %>%  # turn Site into an ordered factor
  ordr::ordinate(
    cols = any_of(spectraList),
    model = ~ prcomp(clr(.)),
    augment = any_of(c("Depth", "Site"))
  ) %>%
  ordr::ggbiplot(., mapping = aes(x = 1, y = 2), sec.axes = "cols", scale.factor = 3 #, groups = CONISS, ellipse = TRUE, ellipse.prob = 0.95, circle = TRUE
  ) +
  ordr::geom_rows_point(aes(colour = Site), alpha = 0.7, size = 2) +
  #ordr::geom_cols_vector() +
  #ordr::geom_cols_text_radiate(aes(label = name), size = 3) +
  scale_color_manual(values = site_colours) + # use custom colour scale
  #   scale_color_viridis_d(option = "C") +  # default colour scale 
  labs(title = "PCA of all cores and universal model (spectra)", color = "Site") +
  theme_minimal() +
  theme(legend.position = "right")


# Create PCA plot with both colour and shape aesthetics
VNIR_Unimodel %>%
  arrange(site_order) %>%
  mutate(Site = factor(Site, levels = unique(Site))) %>%
  ordr::ordinate(
    cols = any_of(spectraList),
    model = ~ prcomp(clr(.)),
    augment = any_of(c("Depth", "Site", "Type"))
  ) %>%
  ordr::ggbiplot(., mapping = aes(x = 1, y = 2), sec.axes = "cols", scale.factor = 3) +
  ordr::geom_rows_point(aes(colour = Site, shape = Type), alpha = 0.7, size = 2) +
  scale_color_manual(values = site_colours) +
  labs(
    title = "PCA of all cores and universal model (spectra)",
    color = "Site",
    shape = "Type") +
  theme_minimal() +
  theme(legend.position = "right")

# PC2 vs PC3
VNIR_Unimodel %>%
  arrange(site_order) %>%  # Order rows by custom numeric order
  mutate(Site = factor(Site, levels = unique(Site))) %>%  # Ordered factor for Site
  ordr::ordinate(
    cols = any_of(spectraList),
    model = ~ prcomp(clr(.)),
    augment = any_of(c("Depth", "Site", "Type"))
  ) %>%
  ordr::ggbiplot(mapping = aes(x = 2, y = 3), sec.axes = "cols", scale.factor = 3) +  # PC2 vs PC3
  ordr::geom_rows_point(aes(colour = Site, shape = Type), alpha = 0.7, size = 2) +
  scale_color_manual(values = site_colours) +
  labs(title = "PCA of Spectra (PC2 vs PC3)",
       color = "Site",
       shape = "Type") +
  theme_minimal() +
  theme(legend.position = "right")





### 5.4 Diagenetic layer ----

# As mentioned in section 3.2, the top layer in sediments is affected by diagenesis.
# Here we perform a PCA to show the diagenetic-sensitive elements driving the variability.

XRF_avg %>%
  mutate(layer = ifelse(Depth <= 3, 
                        "Top layer (0–3 cm)", 
                        "Below 3 cm")) %>%   # create depth layer column
  arrange(site_order) %>%  
  mutate(Site = factor(Site, levels = unique(Site))) %>%  
  ordr::ordinate(
    cols = any_of(elementsList),
    model = ~ prcomp(clr(.)),
    augment = any_of(c("Depth", "Site", "layer"))
  ) %>%
  ordr::ggbiplot(sec.axes = "cols", scale.factor = 8) +
  ordr::geom_rows_point(aes(colour = layer), alpha = 0.7, size = 2) +
  ordr::geom_cols_vector() +
  ordr::geom_cols_text_radiate(aes(label = name), size = 3) +
  scale_color_manual(values = c("Top layer (0–3 cm)" = "deeppink3", 
                                "Below 3 cm" = "darkkhaki")) +
  labs(title = "PCA of all elements for all sites", color = "Depth layer") +
  theme_minimal() +
  theme(legend.position = "right")





### 5.5 Elemental pie chart ----

# to get an overview of average sediment composition

XRF_avg %>%
  summarise(across(5:29, mean, na.rm = TRUE)) %>%
  pivot_longer(cols = everything(), names_to = "Element", values_to = "Mean") %>%
  arrange(desc(Mean)) %>%
  mutate(
    Element = factor(Element, levels = Element),
    Perc = Mean / sum(Mean) * 100
  ) %>%
  ggplot(aes(x = "", y = Mean, fill = Element)) +
  geom_bar(stat = "identity") +
  coord_polar(theta = "y") +
  geom_text(
    aes(label = ifelse(Perc > 2, paste0(round(Perc, 1), "%"), "")), # only show > 2%
    position = position_stack(vjust = 0.5),
    size = 3
  ) +
  scale_fill_manual(values = colorRampPalette(brewer.pal(12, "Set3"))(25)) + # auto 25 colours
  labs(
    title = "Average Elemental Composition (XRF)",
    fill = "Element"
  ) +
  theme_void() +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    legend.title = element_text(face = "bold")
  )





### 5.6 Example code for saving plots (as pdf) ----


# The following scripts can be used to save desired plots in pdf format.
# First create the pdf, add title, execute plot, then save/close the pdf.
# When using models, it can be a good idea to print the summary into the pdf. 
# Instructions down below show how to add the summary to pdf in 2 different ways, 
# depending on the length of the summary.


# # # #

pdf(file = "~/Documents/Sorbonne/Master/Stage M2/data/SedimentData/VNIR_pca.pdf", width = 9, height = 6)
# save as pdf file

text1=paste("PCA for VNIR spectra 
            for all cores in the Gulf of Bothnia")
### positioning of text explaining
plot(NA, xlim=c(0,5), ylim=c(0,5), bty='n', xaxt='n', yaxt='n', xlab='', ylab='')
text(1,4,text1, pos=4)

# # # #


# execute chosen plot


# # # #

dev.off() # closes/saves pdf. To use at very end.

# # # #




# if you want to print the summary of a model in the pdf, use the following:

# # # # # # # # # # # # #

pdf(file = "~/Documents/Sorbonne/Master/Stage M2/data/TimeSeries/ChlaTS.pdf", width = 9, height = 6)
# save as pdf file

text1=paste("Time series for Chla")
# positioning of text explaining
plot(NA, xlim=c(0,5), ylim=c(0,5), bty='n', xaxt='n', yaxt='n', xlab='', ylab='')
text(1,4,text1, pos=4)

# # # # # # # # # # # # #


# execute chosen plot


# New page for summary
plot.new()  # starts a new blank page
text(0, 1, "Model Summary:", pos = 4, cex = 1.2)


# Capture summary and print line by line
summary_lines <- capture.output(summary(MODELNAME))

# Add the summary to the PDF
for (i in seq_along(summary_lines)) {
  text(0, 0.95 - i * 0.035, summary_lines[i], adj = 0, cex = 0.7)
}


# if summary too long, use following code:

summary_lines <- capture.output(summary(MODELNAME))
lines_per_page <- 25  # adjust depending on font size
num_pages <- ceiling(length(summary_lines) / lines_per_page)

for (p in seq_len(num_pages)) {
  plot.new()
  text(0, 1, paste("Model Summary - Page", p), pos = 4, cex = 1.2)
  
  start <- (p - 1) * lines_per_page + 1
  end <- min(p * lines_per_page, length(summary_lines))
  
  for (i in seq_along(summary_lines[start:end])) {
    text(0, 0.95 - i * 0.035, summary_lines[start:end][i], adj = 0, cex = 0.7)
  }
}


# # # #

dev.off() # closes/saves pdf. To use at very end.

# # # #





# # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # 
#
# End of script
#
# # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # 