#Casorso et al. 2026 - R Code Associated with Analyses

#This code is broken up into six overarching segments and is intended to be run in order. To collapse all segments, use Alt + O (Windows/Linux) or Option + O (macOS).
#If you want to run one analysis section independently, you will need to run Code Segments 0 and 1 first. (Note Sections 3.1 to 3.2 and Sections 4.1 to 4.4 each must be run together in order.)
#Intended for use in R version 4.6.0 (or later).


#Code Segment 0) Preparation--------------

#Clear your environment
rm(list=ls())

#Install/Load required packages
  #List all required packages:
packages <- c("dplyr","stringr","ape","caper","tidyr","sf","MODISTools","MuMIn","glmmTMB","ggplot2","readr","car","betareg","AICcmodavg","emmeans","effectsize")
  #Install any missing packages and load all:
package_load <- lapply(packages, function(pkg) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
})

  #Install/Load V.PhyloMaker2 package from GitHub repository
#install.packages("devtools")
#library(devtools)
#devtools::install_github("jinyizju/V.PhyloMaker2")
library(V.PhyloMaker2)


###
#Code Segment 1) Dataset filtering------

#Load raw dataset from working directory ('fileEncoding = "latin1"' prevents encoding errors)
setwd("C:/Users/Casorso/Desktop")
df<-read.csv("Casorso_et_al_2026_-_Raw_Dataset.csv", fileEncoding = "latin1")

#Assign relevant variables as factors or numeric
df$RipenessStage <- as.factor(df$RipenessStage)
df$SampleType <- as.factor(df$SampleType)
df$ExcludeFromFinalAnalysis <- as.factor(df$ExcludeFromFinalAnalysis)
df$SpeciesMin3MRRORFruits <- as.factor(df$SpeciesMin3MRRORFruits)
df$FruitSpeciesName <- as.factor(df$FruitSpeciesName)
df$FruitUniqueID <- as.factor(df$FruitUniqueID)
df$SampleUniqueID <- as.factor(df$SampleUniqueID)
df$UniqueMeasurementID <- as.factor(df$UniqueMeasurementID)
df$ProportionEtOH <- as.numeric(df$ProportionEtOH, na.rm=TRUE)
df$BreathalyzerMeasurement <- as.numeric(df$BreathalyzerMeasurement, na.rm=TRUE)

#Create new columns for unique spatial points (i.e. unique latitude/longitude coordinate sets), and for unique fruit species + spatial points (i.e. unique combinations of fruit species and lat/long coordinate sets); these will be used when controlling for sampling location in later analyses 
df <- df %>%
  mutate(PointID = paste0("Point_", Latitude, "_", Longitude)) %>% #unique spatial points
  mutate(FruitPointID = paste0(FruitSpeciesName, "_", PointID)) #unique fruit species + spatial points


  #Filter to obtain the Full Fruit Dataset 

#Remove unreliable records, species not classified to at least taxonomic genus, and floral tissue samples:
df_full_fruit <- df %>% 
  filter(ExcludeFromFinalAnalysis=="No" | ExcludeFromFinalAnalysis=="Presence/absence only") %>% #include only records not to be excluded and "presence/absence only" records (i.e. records without a reliable %-EtOH value that were retained solely as evidence of ethanol presence)
  filter(SampleType=="fruit") %>% #include only fruits
  filter(!grepl("unknown", FruitSpeciesName)) %>% #remove unknown species
  filter(!str_detect(FruitSpeciesName, "aceae")) #remove fruits classified only to taxonomic family

#Save dataset to working directory, if desired:
#write.csv(df_full_fruit, "Full_fruit_dataset.csv")


##


  #Filter to obtain Repeatedly Sampled Dataset

#Filter according to different minimum sampling thresholds, AND exclude unripe fruits and other ripeness stages accordingly:
#1) include only species with at least three ripe (R) fruits sampled (i.e. '3R')
df_3R <- df_full_fruit %>%
  filter(SpeciesMin3RFruits=="Yes") %>% #include only species that meet sampling threshold
  filter(RipenessStage=="ripe") #include only ripe fruits

#2) include only species with at least three ripe (R) or overripe (OR) fruits sampled (i.e. '3ROR')
df_3ROR <- df_full_fruit %>%
  filter(SpeciesMin3RORFruits=="Yes") %>% #include only species that meet sampling threshold
  filter(RipenessStage=="ripe" | RipenessStage=="overripe") #include only ripe and overripe fruits

#3) include only species with at least three midripe (MR), ripe (R), or overripe (OR) fruits sampled (i.e. '3MRROR')
df_3MRROR <- df_full_fruit %>%
  filter(SpeciesMin3MRRORFruits=="Yes") %>% #include only species that meet sampling threshold
  filter(RipenessStage=="midripe" | RipenessStage=="ripe" | RipenessStage=="overripe") #include only midripe, ripe, and overripe fruits

#Save datasets to working directory folder, if desired:
#write.csv(df_3R, "Repeatedly_sampled_dataset_3-R.csv")
#write.csv(df_3ROR, "Repeatedly_sampled_dataset_3-R-OR.csv")
#write.csv(df_3MRROR, "Repeatedly_sampled_dataset_3-MR-R-OR.csv")


#We conducted Phylogenetic, Environmental, and Canopy vs Ground analyses separately on each of the three sampling thresholds. Results were qualitatively similar across the thresholds, therefore in the manuscript we only focus on and report results for the least conservative threshold (i.e. 'df_3MRROR') and treat this as our "Repeatedly Sampled Dataset". Note that any numerical results included in code notes in this R script correspond to this threshold, as well (excluding those code chunks corresponding to the full fruit dataset).
#At the start of each of these analyses' code segments, optional lines of code are included to instead use one of the other stricter thresholds (i.e. 'df_3R' or 'df_3ROR') as the 'repeatedly sampled dataset' so that comparative analyses can be carried out if desired. 


###
#Code Segment 2) Summary Statistics of Fruits & Floral Tissues--------

#This segment allows calculation of the various descriptive summary statistics reported in the manuscript (including descriptive analysis of floral tissue samples)


#1) Ethanol trends across all fruits (i.e. full fruit dataset ("FF"); all species regardless of sampling repetition and ripeness stages)

  #1a) Number of fruit species with ethanol detected:
#Filter to obtain subset of only the fruit species that had ethanol detected
df_etoh_FF <- df_full_fruit %>%
  filter(EthanolDetectedInSpecies == "Yes") #include only species with ethanol detected

#Calculate percentage of species with ethanol detected
n_EtOH_detected_FF <- length(unique(df_etoh_FF$FruitSpeciesName)) #number of species with ethanol detected (= 272)
n_all_fruits_FF <- length(unique(df_full_fruit$FruitSpeciesName)) #total number of fruit species (in full fruit dataset) (= 591)
n_EtOH_detected_FF / n_all_fruits_FF *100
#ethanol detected in 46.0% of all fruit species 


  #1b)  Measures of central tendency
df_full_fruit %>%
  filter(ProportionEtOH != 0) %>% #remove measurements of 0% ethanol
  summarize(lower_bound = quantile(ProportionEtOH*100, 0.05, na.rm = TRUE), #calculate 5th percentile
            upper_bound = quantile(ProportionEtOH*100, 0.95, na.rm = TRUE), #calculate 95th percentile
            mean = mean(ProportionEtOH*100, na.rm=TRUE), #calculate overall mean value
            median = median(ProportionEtOH*100, na.rm=TRUE)) #calculate overall median value
#the central 90% of values fall between 0.03% and 1.05%
#mean = 0.26%; median = 0.12%


  #1c) Ethanol distributions across ripeness stages

#Calculate maximum ethanol concentrations detected for each ripeness stage:
aggregate(x=df_full_fruit$ProportionEtOH*100, by=list(df_full_fruit$RipenessStage), na.rm=TRUE, FUN=max)
#unripe = 1.49%; midripe = 3.00%; ripe = 3.75%; overripe = 3.90%

#Plot ethanol distributions across ripeness stages in species with ethanol detected (Supplementary Figure S2):
df_ripeness_plot <- df_etoh_FF %>%
  filter(ExcludeFromFinalAnalysis=="No" ) %>% #include only records with reliable %-EtOH values
  filter(RipenessStage=="unripe" | RipenessStage=="midripe" | RipenessStage=="ripe" | RipenessStage=="overripe") #this excludes any 'unknown' records
ripeness_order <- factor(df_ripeness_plot$RipenessStage, level = c("unripe", "midripe", "ripe", "overripe")) #re-order ripeness stage for the plot
plot_ripeness <- ggplot(data = df_ripeness_plot, aes(x = ripeness_order, y = ProportionEtOH*100)) + 
  geom_boxplot(outlier.shape=NA) + 
  geom_jitter(position = position_jitter(w=0.02, h=0), shape=1) + 
  xlab("Ripeness Stage") + ylab("Ethanol Concentration (%)") + theme_bw() ; plot_ripeness 


  #1d) Top ten species with the highest ethanol concentrations detected (Table 1)
top_10_species <- df_full_fruit %>%
  group_by(FruitSpeciesName) %>%
  slice_max(order_by = ProportionEtOH, n = 1, with_ties = FALSE) %>% #find the single highest EtOH value per species
  ungroup() %>%
  arrange(desc(ProportionEtOH)) %>% slice_head(n = 10) %>% #arrange in descending order and select the top 10
  dplyr::select(FruitSpeciesName, ProportionEtOH, RipenessStage) #select relevant columns for final output
top_10_species


##


#2) Ethanol trends across repeatedly-sampled fruits (i.e. repeatedly-sampled dataset ("RS"); species-level sampling repetition and excluding unripe fruits)

#First select a sampling threshold (obtained in Code Segement 1); results reported here and in the manuscript are based on Option 3 
  #Option 1) species with at least 3 ripe fruits: 
    #run code:     
      #df_repeated_samples <- df_3R
  #Option 2) species with at least 3 ripe or overripe fruits:
    #run code:     
      #df_repeated_samples <- df_3ROR
  #Option 3) species with at least 3 midripe, ripe, or overripe fruits: 
    df_repeated_samples <- df_3MRROR


  #2a) Number of fruit species with ethanol detected:

#Filter to obtain only species with ethanol detected
df_etoh_RS <- df_repeated_samples %>%
  filter(EthanolDetectedInSpecies == "Yes") #include only species with ethanol detected

#Calculate percentage of species with ethanol detected
n_EtOH_detected_RS <- length(unique(df_etoh_RS$FruitSpeciesName)) #number of species with ethanol detected (= 151)
n_all_fruits_RS <- length(unique(df_repeated_samples$FruitSpeciesName)) #total number of fruit species (in repeatedly sampled dataset) (= 211)
n_EtOH_detected_RS / n_all_fruits_RS *100
#ethanol detected in 71.6% of repeatedly-sampled fruit species   


  #2b)  Measures of central tendency
df_repeated_samples %>%
  filter(ProportionEtOH != 0) %>% #remove measurements of 0% ethanol
  summarize(lower_bound = quantile(ProportionEtOH*100, 0.05, na.rm = TRUE), #calculate 5th percentile
            upper_bound = quantile(ProportionEtOH*100, 0.95, na.rm = TRUE), #calculate 95th percentile
            mean = mean(ProportionEtOH*100, na.rm=TRUE), #calculate overall mean value
            median = median(ProportionEtOH*100, na.rm=TRUE)) #calculate overall median value
#the central 90% of values fall between 0.03% and 1.09%
#mean = 0.27%; median = 0.12%


  #2c) Ethanol distributions across ripeness stages

#Calculate maximum ethanol concentrations detected for each ripeness stage:
aggregate(x=df_repeated_samples$ProportionEtOH*100, by=list(df_repeated_samples$RipenessStage), na.rm=TRUE, FUN=max)
#exactly the same as the full fruit dataset

#Calculate ethanol detection per ripeness stage
detection_ripeness <- df_repeated_samples %>%
  group_by(RipenessStage) %>%
  summarise(total_fruits = length(unique(FruitUniqueID)), #total number unique fruits
            non_zero_fruits = length(unique(FruitUniqueID[BreathalyzerMeasurement != 0])), #number unique fruits with at least one breathalyzer measurement > 0
            percentage = non_zero_fruits / total_fruits *100)
detection_ripeness
#ethanol detected in 33.6% of midripe fruits, 47.0% of ripe fruits, and 69.7% of overripe fruits


##


#3) Floral tissue samples

#Filter to obtain subset of only floral tissue samples
df_floral <- df %>%
  filter(ExcludeFromFinalAnalysis=="No" | ExcludeFromFinalAnalysis=="Presence/absence only") %>% #include only records not to be excluded and "presence/absence only" records (i.e. records without a reliable %-EtOH value that were retained solely as evidence of ethanol presence)
  filter(SampleType=="flower" | SampleType=="flower & fruit" | SampleType=="nectar" | SampleType=="bract") #include floral tissue samples 


  #3a) Ethanol trends in flowers (mature & flower buds) only

#Filter for only the flower samples
df_flowers <- df_floral %>%
  filter(SampleSubtype != "liquid in bract") %>% filter(SampleSubtype != "nectar") #remove non-flower samples

#Filter for flower samples with ethanol detected
df_flowers_etoh <- df_flowers %>%
  filter(!is.na(BreathalyzerMeasurement), BreathalyzerMeasurement > 0)

#Calculate percentage of species with ethanol detected
n_EtOH_detected_flower <- length(unique(df_flowers_etoh$FruitSpeciesName)) #number of unique species with ethanol detected (= 12)
n_all_flower <- length(unique(df_flowers$FruitSpeciesName)) #total number of floral species (in full fruit dataset) (= 37)
n_EtOH_detected_flower / n_all_flower *100
#ethanol detected in 32.4% of all fruit species 

#Calculate minimum and maximum non-zero ethanol concentrations
min(df_flowers_etoh$ProportionEtOH*100, na.rm=TRUE) #minimum = 0.04%
max(df_flowers_etoh$ProportionEtOH*100, na.rm=TRUE) #maximum = 0.37%


  #3b) Ethanol trends in floral exudates (nectar and liquid in bracts) only

#Filter for only the floral exudate samples (only 1 species)
df_exudates <- df_floral %>%
  filter(SampleSubtype=="liquid in bract" | SampleSubtype=="nectar") #include only exudate samples

#Filter for exudate samples with ethanol detected
df_exudates_etoh <- df_exudates %>%
  filter(!is.na(BreathalyzerMeasurement), BreathalyzerMeasurement > 0)

#Calculate maximum non-zero ethanol concentrations in each sample type
max(df_exudates_etoh$ProportionEtOH[df_exudates_etoh$SampleSubtype=="liquid in bract"]*100, na.rm=TRUE) #maximum for bract liquid = 0.43%
max(df_exudates_etoh$ProportionEtOH[df_exudates_etoh$SampleSubtype=="nectar"]*100, na.rm=TRUE) #maximum for nectar = 0.73%


  #3c) Plot ethanol distribution across all floral tissue sample types (Supplementary Figure S3):
df_floral_plot <- df_floral %>% 
  filter(ExcludeFromFinalAnalysis=="No") %>% #include only records with reliable %-EtOH values
  #Rename a couple groups for purposes of simplifying the plot:
  mutate(SampleSubtype = recode_values(SampleSubtype,
                                       "mature flowers + fruits" ~ "mature flower", #this sample was primarily mature flowers
                                       "mature flowers + flower buds" ~ "mature flower",
                                       "bracts (immature flowers)" ~ "flower bud",
                                       "mature flowers" ~ "mature flower",
                                       default = SampleSubtype))
plot_floral <- ggplot(data = df_floral_plot, aes(x = reorder(SampleSubtype,ProportionEtOH*100,FUN=median), y = ProportionEtOH*100)) +
  geom_boxplot(outlier.shape=NA) + 
  geom_jitter(position = position_jitter(w=0.02, h=0), shape=1) + 
  xlab("Sample Type") + ylab("Ethanol Concentration (%)") + theme_bw() ; plot_floral


###
#Code Segment 3.1) Phylogenetic Tree Creation--------------

#Step 1: Import species lists (there are unique lists for each sampling threshold, as well as for ALL fruit species regardless of sampling):
sp_list_R <- read.csv("Casorso_et_al_2026_-_Phylo_Tree_Species_List_3R.csv") #species with at least 3 ripe fruits
sp_list_ROR <- read.csv("Casorso_et_al_2026_-_Phylo_Tree_Species_List_3ROR.csv") #species with at least 3 ripe or overripe fruits
sp_list_MRROR <- read.csv("Casorso_et_al_2026_-_Phylo_Tree_Species_List_3MRROR.csv") #species with at least 3 midripe, ripe, or overripe fruits
sp_list_all <- read.csv("Casorso_et_al_2026_-_Phylo_Tree_Species_List_All.csv") #no sampling threshold

#Choose the sampling threshold to use going forward ('RS' = Repeatedly Sampled dataset):
  #(for our final analysis, we used "sp_list_MRROR")
#sp_list_RS <- sp_list_R
#sp_list_RS <- sp_list_ROR
sp_list_RS <- sp_list_MRROR

#Repeat for all species, with no sampling threshold and all ripeness stages included ('FF' = Full Fruit Dataset):
sp_list_FF <- sp_list_all


#Step 2: Generate phylogenetic trees
  
  #(a) Repeated Sampling:
tree_output_RS <- phylo.maker(sp.list = sp_list_RS, tree = GBOTB.extended.WP, nodes=nodes.info.1.WP) #typically takes a minute to run
tree_RS <- tree_output_RS$scenario.3 #extract tree
  
  #(b) All species:
tree_output_FF <- phylo.maker(sp.list = sp_list_FF, tree = GBOTB.extended.WP, nodes=nodes.info.1.WP) #typically takes a minute to run
tree_FF <- tree_output_FF$scenario.3 #extract tree

#Save tree newick strings to working directory folder, if desired (these can be used to visualize trees by uploading to a tree visualization program, such as iTOL, which was used to create the figures in the manuscript)
#write.tree(tree_RS, "PhylogeneticTree_RepeatedlySampledDataset.nwk")
#write.tree(tree_FF, "PhylogeneticTree_FullFruitDataset.nwk)


###
#Code Segment 3.2) Phylogenetic Analysis (Pagel's Lambda)--------------

#First we will carry out the analysis using the repeatedly sampled dataset (or an alternative minimum sampling threshold), then repeat for the full fruit dataset (i.e. all species & ripeness stages)


  #Step 1) Calculate species metrics (median, mean, and maximum EtOH)

#First select the appropriate sampling threshold (these are obtained from Code Segment 1); for final analysis, we used Option 3 
  #Option 1) species with at least 3 ripe fruits: 
    #run code:     
      #df_repeated_samples <- df_3R
  #Option 2) species with at least 3 ripe or overripe fruits:
    #run code:     
      #df_repeated_samples <- df_3ROR
  #Option 3) species with at least 3 midripe, ripe, or overripe fruits: 
    df_repeated_samples <- df_3MRROR

#Calculate metrics:
#Note we need to remove records maintained in the repeatedly sampled dataset which do not have a reliable % ethanol value (i.e. records maintained solely as evidence of ethanol presence/absence)
sp_metrics_phylo <- df_repeated_samples %>%
  filter(ExcludeFromFinalAnalysis=="No") %>% #include only records with reliable %-EtOH values
  group_by(FruitSpeciesName_phylo) %>%
  summarise(MedianEtOH = median(ProportionEtOH), #median ethanol
            MeanEtOH = mean(ProportionEtOH), #mean ethanol
            MaximumEtOH = max(ProportionEtOH)) #maximum ethanol
sp_metrics_phylo <- as.data.frame(sp_metrics_phylo) #write output into a dataframe

    
  #Step 2) Calculate lambda with a maximum likelihood (ML) estimate, for each metric

#First combine phylogeny and metric data into a single object
#Use tree 'tree_RS' from Code Segment 3.1, **BUT** ensure the same sampling threshold was used to create both the tree AND the species-specific metrics.
tree_RS$node.label<-NULL #ignore node labels
EtOHPhyloData <- comparative.data(phy = tree_RS, data = sp_metrics_phylo, names.col = FruitSpeciesName_phylo, vcv = TRUE, na.omit = FALSE, warn.dropped = TRUE)
EtOHPhyloData$dropped #check if any species were dropped (should not have any dropped)

  #(a) Median:
median_lambda <- pgls(MedianEtOH ~ 1, data = EtOHPhyloData, lambda='ML')
summary(median_lambda)
#lambda = 0.000
#lambda is significantly different from 1 (p < 0.05), and is not significantly different from 0 (p > 0.05) --> consistent with no phylogenetic signal

  #(b) Mean:
mean_lambda <- pgls(MeanEtOH ~ 1, data = EtOHPhyloData, lambda='ML')
summary(mean_lambda)
#lambda = 0.000
#lambda is significantly different from 1 (p < 0.05), and is not significantly different from 0 (p > 0.05) --> consistent with no phylogenetic signal

  #(c) Maximum:
maximum_lambda <- pgls(MaximumEtOH ~ 1, data = EtOHPhyloData, lambda='ML')
summary(maximum_lambda)
#lambda = 0.697
#lambda is significantly different from 0 (p < 0.05) AND from 1 (p < 0.05) --> consistent with intermediate phylogenetic signal


##


#Re-run above phylogenetic analysis on full fruit dataset (i.e. all species and ripeness stages, no minimum sampling threshold)
  #We will use the Full Fruit Dataset ('df_full_fruit') and the full phylogenetic tree ('tree_FF')

#Calculate metrics:
sp_metrics_phylo_FF <- df_full_fruit %>%
  filter(ExcludeFromFinalAnalysis=="No") %>% #include only records with reliable %-EtOH values
  group_by(FruitSpeciesName_phylo) %>%
  summarise(MedianEtOH = median(ProportionEtOH, na.rm=TRUE), #median ethanol
            MeanEtOH = mean(ProportionEtOH, na.rm=TRUE), #mean ethanol
            MaximumEtOH = max(ProportionEtOH, na.rm=TRUE)) #maximum ethanol
sp_metrics_phylo_FF <- as.data.frame(sp_metrics_phylo_FF) #write the output into a dataframe


#Combine phylogenetic tree and metric data into a single object
tree_FF$node.label<-NULL #ignore node labels
EtOHPhyloData_FF <- comparative.data(phy = tree_FF, data = sp_metrics_phylo_FF, names.col = FruitSpeciesName_phylo, vcv = TRUE, na.omit = FALSE, warn.dropped = TRUE)
EtOHPhyloData_FF$dropped #check if any species were dropped (should not have any dropped)


#Calculate Pagel's lambda for each metric
  #(a) Median:
median_lambda_FF <- pgls(MedianEtOH ~ 1, data = EtOHPhyloData_FF, lambda='ML')
summary(median_lambda_FF)
#lambda = 0.000
#lambda is significantly different from 1 (p < 0.05), and is not significantly different from 0 (p > 0.05) --> consistent with no phylogenetic signal
  
  #(b) Mean:
mean_lambda_FF <- pgls(MeanEtOH ~ 1, data = EtOHPhyloData_FF, lambda='ML')
summary(mean_lambda_FF)
#lambda = 0.000
#lambda is significantly different from 1 (p < 0.05), and is not significantly different from 0 (p > 0.05) --> consistent with no phylogenetic signal

  #(c) Maximum:
maximum_lambda_FF <- pgls(MaximumEtOH ~ 1, data = EtOHPhyloData_FF, lambda='ML')
summary(maximum_lambda_FF)
#lambda = 0.496
#lambda is significantly different from 0 (p < 0.05) AND from 1 (p < 0.05) --> consistent with intermediate phylogenetic signal


###
#Code Segment 4.1) Extracting Environmental Data-------------

#Note we will create a temporary working directory folder to store the (very large) environmental data, which will be wiped after closing the RStudio session (working directory location was set in Code Segment 1).

#We will extract and bind environmental data to all sampling locations from the Full Fruit Dataset ('df_full_fruit'), and apply the minimum sampling threshold filters afterwards (i.e. repeatedly sampled dataset), therefore avoiding having to extract & bind the environmental data multiple times to different data subsets. 

#First assemble a dataframe of selected attributes for each sampling location ("sites")
data_sites <- distinct(df_full_fruit %>% dplyr::select(PointID, Country, Latitude, Longitude)) %>%
  filter(!is.na(Latitude))
NPoints <- nrow(data_sites)


##


#Extract and bind selected environmental data from the WorldClim database version 2.1

  #load the database:
rast_worldclim <- geodata::worldclim_global(var="bio", res=5, path= tempdir(), version="2.1")

  #add selected variables as new new columns in the dataframe: 
data_sites$bio1 <- terra::extract(rast_worldclim[["wc2.1_5m_bio_1"]], cbind(data_sites$Longitude, data_sites$Latitude))[,1]  ## BIO1 = Annual Mean Temperature (degrees C)
data_sites$bio12 <- terra::extract(rast_worldclim[["wc2.1_5m_bio_12"]], cbind(data_sites$Longitude, data_sites$Latitude))[,1]  ## BIO12 = Annual Precipitation (mm)
data_sites$bio15 <- terra::extract(rast_worldclim[["wc2.1_5m_bio_15"]], cbind(data_sites$Longitude, data_sites$Latitude))[,1]  ## BIO15 = Precipitation Seasonality (Coefficient of Variation)


##


#Extract and bind satellite NDVI data from MODIS MOD13Q1 database version 6

#First create empty columns in the dataframe for the NDVI mean and SD values that we can then fill in with the data from MODIS
data_sites$NDVI.mean <- rep(NA< NPoints)
data_sites$NDVI.SD <- rep(NA< NPoints)

#Next extract the data and fill in the empty columns; this code will run multiple loops if download interrupts (note this code takes several minutes to run)
for(x in 1:(NPoints*2)){
  p <- which(is.na(data_sites$NDVI.mean))[1]
  if(!is.na(p)){
    start.date = "2023-01-01"       
    end.date = "2023-12-30"   
    buffer_p <- st_buffer(sf::st_point(c(data_sites$Longitude[p],data_sites$Latitude[p])), dist=10/111.1)
    # Get MODIS data from MOD13Q1
    tryCatch(point_ndvi <- mt_subset(product = "MOD13Q1",
                                     lat = data_sites$Latitude[p],
                                     lon =  data_sites$Longitude[p],
                                     band = "250m_16_days_NDVI",
                                     start = start.date,
                                     end = end.date,
                                     km_lr = 1,
                                     km_ab = 1,
                                     site_name = "p",
                                     internal = TRUE,
                                     progress = FALSE), finally = print("Success"), error = function(e) { skip_to_next <<- TRUE})
    if(exists("point_ndvi")){
      data_sites$NDVI.mean[p] <- round(mean(point_ndvi$value, na.rm=T))
      data_sites$NDVI.SD[p] <- round(sd(point_ndvi$value))
      print(paste("datpoint:", p, "   loop:", x))
    }else{}
  }else{}
  rm(point_ndvi)
  if(p==NPoints){break}else{}
}

#Apply a scaling factor to the values
data_sites$NDVI.mean <- data_sites$NDVI.mean*0.0001
data_sites$NDVI.SD <- data_sites$NDVI.SD*0.0001


##


#Extract and bind biogeographical realms data

#First, the map of realms/regions must be downloaded from here: <https://macroecology.ku.dk/resources/wallace>
  #Select the 'Download all maps' tab under 'GIS maps of CMEC zoogeographic realms and regions'.
  #Once the .zip file has downloaded, extract all the files to your working directory (the files should be directly in your working directory location - not within another folder - or else you will need to update your working directory accordingly).


#Next, prepare site coordinates to facilitate binding the realms/regions data
poly_Realms <- st_read("newRealms.shp") #read the 'newRealms' vector shape file into R
pts.multi_all <- sf::st_multipoint(as.matrix(data.frame(data_sites$Longitude, data_sites$Latitude))) #convert latitude/longitude coordinates into a coordinate matrix and group together into one object
pts <- st_cast(st_sfc(pts.multi_all), "POINT") #convert the object into a list of points
st_crs(pts) <- st_crs(poly_Realms) #copy the Coordinate Reference System (CRS) to the list of points
poly_Realms$Realm[which(poly_Realms$Realm=="Oceanina")] <- "Oceania" #Correct wrong spelling in realm name "Oceania"
realm_name <- unique(poly_Realms$Realm) #list of unique realm names
nrealm <- length(realm_name) #number of unique realm names


#Now bind the realms/regions data to our data sites 
  #create empty column in the dataframe for the realm to be filled in:
data_sites$realm <- rep(NA, NPoints)
  #fill in the realm column:
for(r in 1:nrealm){
  sf_use_s2(FALSE)
  polyrealm <- poly_Realms %>% filter (Realm==realm_name[r])
  polyrealm <-  polyrealm[1]
  st_intersection(pts, polyrealm)
  data_sites$realm[which(st_intersects(pts, polyrealm, sparse = FALSE))] <- realm_name[r]
}
  #replace non-captured locations:
data_sites$realm[data_sites$Country=="Taiwan" & is.na(data_sites$realm)] <- "Oriental"
data_sites$realm[data_sites$Country=="Borneo, Malaysia" & is.na(data_sites$realm)] <- "Oriental"


##


#Save site-level data as .csv file, if desired
#write.csv(data_sites, "SiteSpecificEnvironmentalData.csv")


###
#Code Segment 4.2) Binding Environmental Data to the Ethanol Measurements---------------

#Step 1) Assemble environmental attributes for each unique combination of species + site:

  #Create a subset of the environmental dataset (removing the Latitude, Longitude, and Country columns to avoid duplication when merging later):
data_sites_sub <- dplyr::select(data_sites, - c("Latitude", "Longitude", "Country"))

  #Create a subset of the ethanol dataset: 
df_full_fruit_sub <- df_full_fruit  %>%
  dplyr::select(FruitPointID, PointID, Latitude, Longitude, Country, StudySiteName, FruitSpeciesName) %>%
  group_by(FruitPointID) %>%
  slice(1) %>% #select only 1 row of each FruitPointID
  ungroup()

  #Merge both subsets together: 
site_species_attr <- left_join(df_full_fruit_sub, data_sites_sub, by= "PointID")


#Step 2) Select the sampling threshold to use going forward as the Repeatedly Sampled Dataset (these are obtained from Code Segment 1); for our final analyses, we used Option 3.
  #Option 1) species with at least 3 ripe fruits: 
    #run code:     
      #df_repeated_samples <- df_3R
  #Option 2) species with at least 3 ripe or overripe fruits:
    #run code:     
      #df_repeated_samples <- df_3ROR
  #Option 3) species with at least 3 midripe, ripe, or overripe fruits: 
    df_repeated_samples <- df_3MRROR


#Step 3) Next, calculate plant-species- and sampling-location-specific ethanol metrics, and merge with the environmental data

##Note we need to do this separately for (1) Ethanol Occurrence analysis and (2) Mean/Maximum Ethanol Concentration analysis. For (1) Ethanol Occurrence, we can include records maintained in the repeatedly sampled dataset which do not have a reliable % ethanol value and were retained solely as evidence of ethanol presence/absence. For (2) Mean/Maximum Ethanol Concentration, these presence/absence only values must be removed.

  #(a) Ethanol Occurrence (note can safely ignore warning messages from below code chunk)
df_enviro_metrics_Occur <- df_repeated_samples %>%
  group_by(FruitPointID) %>%
  summarise(EtOH_Occur = any(EthanolDetectedInSample), #if ethanol detected in sample, then this = TRUE
            n_samp = length(EthanolDetectedInSample), #number of samples
            n_positive = length(which(EthanolDetectedInSample=='TRUE'))/length(EthanolDetectedInSample), #proportion of samples with ethanol detected
            p_detect = n_positive/n_samp) #probability of detecting ethanol

#Merge species-site attributes with the calculated metrics
df_enviro_metrics_Occur <- left_join(df_enviro_metrics_Occur, site_species_attr, by="FruitPointID")

#Filter to keep only the first unique row for each distinct fruit species:
df_enviro_metrics_Occur <- df_enviro_metrics_Occur %>% distinct(FruitSpeciesName, .keep_all = TRUE)

#Create a position variable for the autocorrelation structure:   
df_enviro_metrics_Occur$pos <- numFactor(round(df_enviro_metrics_Occur$Longitude), round(df_enviro_metrics_Occur$Latitude))


  #(b) Mean/Maximum Ethanol Concentration ('Conc' = concentration)
df_enviro_metrics_Conc <- df_repeated_samples %>% 
  filter(ExcludeFromFinalAnalysis=="No") %>% #filter to remove records without reliable %-EtOH values
  group_by(FruitPointID) %>%
  summarise(EtOH_mean = mean(ProportionEtOH, na.rm = TRUE), #mean EtOH
            EtOH_max = max(ProportionEtOH, na.rm = TRUE), #maximum EtOH
            n_samp = length(ProportionEtOH), #number of samples
            n_positive = length(which(ProportionEtOH>0))/n_samp, #proportion of samples with non-zero ethanol
            p_detect = n_positive/n_samp) #probability of detecting non-zero ethanol

#Merge species-site attributes with the calculated metrics
df_enviro_metrics_Conc <- left_join(df_enviro_metrics_Conc, site_species_attr, by="FruitPointID")

#Filter to keep only the first unique row for each distinct fruit species:
df_enviro_metrics_Conc <- df_enviro_metrics_Conc %>% distinct(FruitSpeciesName, .keep_all = TRUE)

#Create a position variable for the autocorrelation structure     
df_enviro_metrics_Conc$pos <- numFactor(round(df_enviro_metrics_Conc$Longitude), round(df_enviro_metrics_Conc$Latitude))


##

#Save .csv files of fruit-sample-level data to working directory, if desired 
#write.csv(df_enviro_metrics_Occur, "FruitSampleEnviromnentalData_Occurrence.csv")
#write.csv(df_enviro_metrics_Conc, "FruitSampleEnviromnentalData_concentration.csv")


####
#Code Segment 4.3) Environmental Analysis - Ethanol Occurrence-------------------

##Note: To run this analysis with a different minimum sampling threshold, return to Code Segment 4.2 Step 2, choose a different threshold and re-run Step 3 before returning here. While the results are reasonably consistent across sampling thresholds, it is worth noting that each threshold yields different "best" models depending on the combination of predictor variables considered, and we prioritized focusing on a single threshold and maximizing sample size. The code included below to calculate model predictions corresponds only to the most supported model based on this threshold, matching those reported in the manuscript.

#We will build several Generalized Linear Mixed Models with our variables of interest, modeled as various combinations of linear and second-order polynomials, and then compare these models with corrected Akaike Information Criterion (AICc).
#We will fit a binomial distribution to each model, given the binomial dependent variable (EtOH Detected = true or false).
##Note the term ‘exp(pos + 0|1)’ models the spatial autocorrelation in the data, whereby observations closer together in space are more similar to each other and the correlation between observations decays exponentially as the distance between them increases.

#Run models
glmm_enviro_Occur_1 <- glmmTMB(EtOH_Occur ~ -1 + realm + scale(bio1) + scale(bio12) + scale(bio15) + scale(NDVI.mean) +
                            exp(pos + 0|1), data=df_enviro_metrics_Occur, family=binomial(link = "logit"))

glmm_enviro_Occur_2 <- glmmTMB(EtOH_Occur ~ -1 +  realm + poly(bio1,2) +  poly(bio12,2) +  poly(bio15,2) + poly(NDVI.mean,2) +  
                            exp(pos + 0|1), data=df_enviro_metrics_Occur, family=binomial(link = "logit"))

glmm_enviro_Occur_3 <- glmmTMB(EtOH_Occur ~ 1 + scale(bio1) + scale(bio12) + scale(bio15) + scale(NDVI.mean) +
                            exp(pos + 0|1), data=df_enviro_metrics_Occur, family=binomial(link = "logit"))

glmm_enviro_Occur_4 <- glmmTMB(EtOH_Occur ~ 1 + poly(bio1,2) +  poly(bio12,2) +  poly(bio15,2) + poly(NDVI.mean,2) +  
                            exp(pos + 0|1), data=df_enviro_metrics_Occur, family=binomial(link = "logit"))

glmm_enviro_Occur_5 <- glmmTMB(EtOH_Occur ~ 1 + scale(bio1) + scale(bio15) +
                            exp(pos + 0|1), data=df_enviro_metrics_Occur, family=binomial(link = "logit"))

glmm_enviro_Occur_6 <- glmmTMB(EtOH_Occur ~ 1 + realm + scale(bio1) + scale(bio12) +
                            exp(pos + 0|1), data=df_enviro_metrics_Occur, family=binomial(link = "logit"))

glmm_enviro_Occur_7 <- glmmTMB(EtOH_Occur ~ - 1 + realm + scale(bio1) + 
                            exp(pos + 0|1), data=df_enviro_metrics_Occur, family=binomial(link = "logit"))

glmm_enviro_Occur_8 <- glmmTMB(EtOH_Occur ~ poly(bio1,2) + poly(bio12,2) + poly(NDVI.mean,2) +
                            exp(pos + 0|1), data=df_enviro_metrics_Occur, family=binomial(link = "logit"))

glmm_enviro_Occur_9 <- glmmTMB(EtOH_Occur ~ realm + poly(bio1,2) +  poly(bio12,2) + poly(NDVI.mean,2) +
                            exp(pos + 0|1), data=df_enviro_metrics_Occur, family=binomial(link = "logit"))

glmm_enviro_Occur_10 <- glmmTMB(EtOH_Occur ~ - 1 + realm +  poly(bio1,2) + 
                             exp(pos + 0|1), data=df_enviro_metrics_Occur, family=binomial(link = "logit"))

glmm_enviro_Occur_11 <- glmmTMB(EtOH_Occur ~ - 1 + realm  + 
                             exp(pos + 0|1), data=df_enviro_metrics_Occur, family=binomial(link = "logit"))

glmm_enviro_Occur_12 <- glmmTMB(EtOH_Occur ~ poly(bio1,2) +
                             exp(pos + 0|1) + (realm|1), data=df_enviro_metrics_Occur, family=binomial(link = "logit"))

#Calculate AICc values for each model
AICs_Occur <- MuMIn::AICc(glmm_enviro_Occur_1,glmm_enviro_Occur_2,glmm_enviro_Occur_3,glmm_enviro_Occur_4,glmm_enviro_Occur_5,glmm_enviro_Occur_6,glmm_enviro_Occur_7,glmm_enviro_Occur_8,glmm_enviro_Occur_9,glmm_enviro_Occur_10,glmm_enviro_Occur_11,glmm_enviro_Occur_12)
AICs_Occur[order(AICs_Occur$AIC), ] #order AICcs from lowest to highest
#Model #4 is the most supported (i.e. lowest AIC)

#Print results for the best model
summary(glmm_enviro_Occur_4)
#BIO1, BIO12, BIO15, and NDVI have significant impact on ethanol detection (p < 0.05)


##


#Calculate model predictions (for Model #4)

#Assign alpha and z-scores so we can calculate confidence intervals on the logit scale
alpha <- 0.05
z_score <- qnorm(1 - alpha / 2)  # 95% confidence interval


  #1) BIO1 (Annual Mean Temperature)

#First create dataframe of the predictor variables, filled with raw BIO1 values and mean values of all other predictors included in the model
pred_Occur_bio1 <- data.frame(
  bio1 = df_enviro_metrics_Occur$bio1,                                    
  bio12 = mean(df_enviro_metrics_Occur$bio12, na.rm=T),
  bio15 = mean(df_enviro_metrics_Occur$bio15, na.rm=T),
  NDVI.mean = mean(df_enviro_metrics_Occur$NDVI.mean, na.rm=T),
  pos =  df_enviro_metrics_Occur[,"pos"])

#Calculate predictions for partial response BIO1
predictions_Occur_bio1 <- predict(glmm_enviro_Occur_4, pred_Occur_bio1, type = "link", se.fit = TRUE)
#Lower and upper bounds of the logit scale predictions
pred_Occur_bio1$fit_Occur_lower <- predictions_Occur_bio1$fit - z_score * predictions_Occur_bio1$se.fit
pred_Occur_bio1$fit_Occur_upper <- predictions_Occur_bio1$fit + z_score * predictions_Occur_bio1$se.fit
#Convert logit scale predictions back to probabilities (inverse logit)
pred_Occur_bio1$pred_probOccur_bio1 <- plogis(predictions_Occur_bio1$fit)
pred_Occur_bio1$pred_probOccur_bio1_lower <- plogis(pred_Occur_bio1$fit_Occur_lower )
pred_Occur_bio1$pred_probOccur_bio1_upper <- plogis(pred_Occur_bio1$fit_Occur_upper)

#Plot the predictions with 95% confidence intervals (Figure 3A)
plot_pred_Occur_bio1 <- ggplot(pred_Occur_bio1, aes(x = bio1, y = pred_probOccur_bio1*100)) + 
  #Predicted fit:
  geom_line(color = "#be0032", lwd=2) + #predicted probabilities
  #Confidence intervals:
  geom_ribbon(aes(ymin = pred_probOccur_bio1_lower*100, ymax = pred_probOccur_bio1_upper*100), fill = "grey", alpha = 0.2) + 
  #Further aesthetics:
  labs(x = "Annual Mean Temperature (°C)", 
       y = "Ethanol Presence (% of species sampled)") +
  theme_bw() +
  #Add layer with raw data points:
  geom_point(aes(x=df_enviro_metrics_Occur$bio1,
                 #Shift y-values equal to 0 by +5 and y-values not equal to 0 (i.e. = 100) by -5, and then add vertical jitter (= 5 units), to prevent points overlapping at 0% and 100%:
                 y=ifelse(df_enviro_metrics_Occur$EtOH_Occur*100 == 0, 
                          df_enviro_metrics_Occur$EtOH_Occur*100 + 5, 
                          df_enviro_metrics_Occur$EtOH_Occur*100 - 5)), 
             colour = "steelblue", alpha=0.3, size = 2, position = position_jitter(height=5, width=0)) +
  ylim(0, 100) ; plot_pred_Occur_bio1

#Finding model-predicted maximum values and confidence intervals
  #First name the variables:
x<-pred_Occur_bio1$bio1 #x (bio1)
y<-pred_Occur_bio1$pred_probOccur_bio1 #y (% occurrence)
upperCI<-pred_Occur_bio1$pred_probOccur_bio1_upper #upper confidence interval
lowerCI<-pred_Occur_bio1$pred_probOccur_bio1_lower #lower confidence interval
  #Calculate values:
max_index <- which.max(y) #index at maximum y
max_x <- x[max_index] ; max_x #value of x at max y (= 25 degrees C)
max_y <- y[max_index]*100 ; max_y #max y (= 83.2% occurrence)
upperCI[max_index]*100 #upper CI (= 91.4%)
lowerCI[max_index]*100 #lower CI (= 70.0%)


  #2) BIO12 (Annual Precipitation)

#Create dataframe of the predictor variables, filled with raw BIO12 values and mean values of all other predictors included in the model
pred_Occur_bio12 <- data.frame(
  bio1 = mean(df_enviro_metrics_Occur$bio1, na.rm=T),                                    
  bio12 = df_enviro_metrics_Occur$bio12,
  bio15 = mean(df_enviro_metrics_Occur$bio15, na.rm=T),
  NDVI.mean = mean(df_enviro_metrics_Occur$NDVI.mean, na.rm=T),
  pos =  df_enviro_metrics_Occur[,"pos"])

#Calculate predictions for partial response BIO12
predictions_Occur_bio12 <- predict(glmm_enviro_Occur_4, pred_Occur_bio12, type = "link", se.fit = TRUE)
#Lower and upper bounds of the logit scale predictions
pred_Occur_bio12$fit_Occur_lower <- predictions_Occur_bio12$fit - z_score * predictions_Occur_bio12$se.fit
pred_Occur_bio12$fit_Occur_upper <- predictions_Occur_bio12$fit + z_score * predictions_Occur_bio12$se.fit
#Convert logit scale predictions back to probabilities (inverse logit)
pred_Occur_bio12$pred_probOccur_bio12 <- plogis(predictions_Occur_bio12$fit)
pred_Occur_bio12$pred_probOccur_bio12_lower <- plogis(pred_Occur_bio12$fit_Occur_lower )
pred_Occur_bio12$pred_probOccur_bio12_upper <- plogis(pred_Occur_bio12$fit_Occur_upper)

#Plot the predictions with 95% confidence intervals (Figure 3B)
plot_pred_Occur_bio12 <- ggplot(pred_Occur_bio12, aes(x = bio12, y = pred_probOccur_bio12*100)) + 
  #Predicted fit:
  geom_line(color = "#be0032", lwd=2) + #predicted probabilities
  #Confidence intervals:
  geom_ribbon(aes(ymin = pred_probOccur_bio12_lower*100, ymax = pred_probOccur_bio12_upper*100), fill = "grey", alpha = 0.2) +
  #Further aesthetics:
  labs(x = "Annual Precipitation (mm)", 
       y = "Ethanol Presence (% of species sampled)") +
  theme_bw()  +
  #Add layer with raw data points:
  geom_point(aes(x=df_enviro_metrics_Occur$bio12,
                 #Shift y-values equal to 0 by +5 and y-values not equal to 0 (i.e. = 100) by -5, and then add vertical jitter (= 5 units), to prevent points overlapping at 0% and 100%:
                 y=ifelse(df_enviro_metrics_Occur$EtOH_Occur*100 == 0, 
                          df_enviro_metrics_Occur$EtOH_Occur*100 + 5, 
                          df_enviro_metrics_Occur$EtOH_Occur*100 - 5)), 
             colour = "steelblue", alpha=0.3, size = 2, position = position_jitter(height=5, width=0)) +
  ylim(0, 100) ; plot_pred_Occur_bio12

#Finding model-predicted maximum values and confidence intervals
  #First name the variables (note we are writing over the value names from earlier):
x<-pred_Occur_bio12$bio12 #x (bio12)
y<-pred_Occur_bio12$pred_probOccur_bio12 #y (% occurrence)
upperCI<-pred_Occur_bio12$pred_probOccur_bio12_upper #upper confidence interval
lowerCI<-pred_Occur_bio12$pred_probOccur_bio12_lower #lower confidence interval
  #Calculate values:
max_index <- which.max(y) #index at maximum y
max_x <- x[max_index] ; max_x #value of x at max y (= 3335 mm)
max_y <- y[max_index]*100 ; max_y #max y (= 95.3% occurrence)
upperCI[max_index]*100 #upper CI (= 99.3%)
lowerCI[max_index]*100 #lower CI (= 74.1%)


  #3) BIO15 (Precipitation Seasonality)

#Create dataframe of the predictor variables, filled with raw BIO15 values and mean values of all other predictors included in the model
pred_Occur_bio15 <- data.frame(
  bio1 = mean(df_enviro_metrics_Occur$bio1, na.rm=T),                                    
  bio12 = mean(df_enviro_metrics_Occur$bio12, na.rm=T),
  bio15 = df_enviro_metrics_Occur$bio15,
  NDVI.mean = mean(df_enviro_metrics_Occur$NDVI.mean, na.rm=T),
  pos =  df_enviro_metrics_Occur[,"pos"])

#Calculate predictions for partial response BIO15
predictions_Occur_bio15 <- predict(glmm_enviro_Occur_4, pred_Occur_bio15, type = "link", se.fit = TRUE)
#Lower and upper bounds of the logit scale predictions
pred_Occur_bio15$fit_Occur_lower <- predictions_Occur_bio15$fit - z_score * predictions_Occur_bio15$se.fit
pred_Occur_bio15$fit_Occur_upper <- predictions_Occur_bio15$fit + z_score * predictions_Occur_bio15$se.fit
#Convert logit scale predictions back to probabilities (inverse logit)
pred_Occur_bio15$pred_probOccur_bio15 <- plogis(predictions_Occur_bio15$fit)
pred_Occur_bio15$pred_probOccur_bio15_lower <- plogis(pred_Occur_bio15$fit_Occur_lower )
pred_Occur_bio15$pred_probOccur_bio15_upper <- plogis(pred_Occur_bio15$fit_Occur_upper)

#Plot the predictions with 95% confidence intervals (Figure 3C)
plot_pred_Occur_bio15 <- ggplot(pred_Occur_bio15, aes(x = bio15, y = pred_probOccur_bio15*100)) + 
  #Predicted fit:
  geom_line(color = "#be0032", lwd=2) + #predicted probabilities
  #Confidence intervals:
  geom_ribbon(aes(ymin = pred_probOccur_bio15_lower*100, ymax = pred_probOccur_bio15_upper*100), fill = "grey", alpha = 0.2) + 
  #Further aesthetics:
  labs(x = "Precipitation Seasonality (CV)", 
       y = "Ethanol Presence (% of species sampled)") +
  theme_bw()  +
  #Add layer with raw data points:
  geom_point(aes(x=df_enviro_metrics_Occur$bio15,
                 #Shift y-values equal to 0 by +5 and y-values not equal to 0 (i.e. = 100) by -5, and then add vertical jitter (= 5 units), to prevent points overlapping at 0% and 100%:
                 y=ifelse(df_enviro_metrics_Occur$EtOH_Occur*100 == 0, 
                          df_enviro_metrics_Occur$EtOH_Occur*100 + 5, 
                          df_enviro_metrics_Occur$EtOH_Occur*100 - 5)), 
             colour = "steelblue", alpha=0.3, size = 2, position = position_jitter(height=5, width=0)) +
  ylim(0, 100) ; plot_pred_Occur_bio15

#Finding model-predicted maximum values and confidence intervals
  #First name the variables (note we are writing over the value names from earlier):
x<-pred_Occur_bio15$bio15 #x (bio15)
y<-pred_Occur_bio15$pred_probOccur_bio15 #y (% occurrence)
upperCI<-pred_Occur_bio15$pred_probOccur_bio15_upper #upper confidence interval
lowerCI<-pred_Occur_bio15$pred_probOccur_bio15_lower #lower confidence interval
  #Calculate values:
max_index <- which.max(y) #index at maximum y
max_x <- x[max_index] ; max_x #value of x at max y (= 133.1%)
max_y <- y[max_index]*100 ; max_y #max y (= 88.1% occurrence)
upperCI[max_index]*100 #upper CI (= 98.3%)
lowerCI[max_index]*100 #lower CI (= 49.0%)


  #4) NDVI

#Create dataframe of the predictor variables, filled with raw NDVI mean values and mean values of all other predictors included in the model
pred_Occur_NDVI <- data.frame(
  bio1 = mean(df_enviro_metrics_Occur$bio1, na.rm=T),                                    
  bio12 = mean(df_enviro_metrics_Occur$bio12, na.rm=T),
  bio15 = mean(df_enviro_metrics_Occur$bio15, na.rm=T),
  NDVI.mean = df_enviro_metrics_Occur$NDVI.mean,
  pos =  df_enviro_metrics_Occur[,"pos"])

#Calculate predictions for partial response NDVI
predictions_Occur_NDVI <- predict(glmm_enviro_Occur_4, pred_Occur_NDVI, type = "link", se.fit = TRUE)
#Lower and upper bounds of the logit scale predictions
pred_Occur_NDVI$fit_Occur_lower <- predictions_Occur_NDVI$fit - z_score * predictions_Occur_NDVI$se.fit
pred_Occur_NDVI$fit_Occur_upper <- predictions_Occur_NDVI$fit + z_score * predictions_Occur_NDVI$se.fit
#Convert logit scale predictions back to probabilities (inverse logit)
pred_Occur_NDVI$pred_probOccur_NDVI <- plogis(predictions_Occur_NDVI$fit)
pred_Occur_NDVI$pred_probOccur_NDVI_lower <- plogis(pred_Occur_NDVI$fit_Occur_lower )
pred_Occur_NDVI$pred_probOccur_NDVI_upper <- plogis(pred_Occur_NDVI$fit_Occur_upper)

#Plot the model predictions with 95% confidence intervals (Figure 3D)
plot_pred_Occur_NDVI <- ggplot(pred_Occur_NDVI, aes(x = NDVI.mean, y = pred_probOccur_NDVI*100)) + 
  #Predicted fit:
  geom_line(color = "#be0032", lwd=2) + 
  #Confidence intervals:
  geom_ribbon(aes(ymin = pred_probOccur_NDVI_lower*100, ymax = pred_probOccur_NDVI_upper*100), fill = "grey", alpha = 0.2) +
  #Further aesthetics:
  labs(x = "NDVI", 
       y = "Ethanol Presence (% of species sampled)") +
  theme_bw()  +
  #Add layer with raw data points:
  geom_point(aes(x=df_enviro_metrics_Occur$NDVI.mean,
                 #Shift y-values equal to 0 by +5 and y-values not equal to 0 (i.e. = 100) by -5, and then add vertical jitter (= 5 units), to prevent points overlapping at 0% and 100%:
                 y=ifelse(df_enviro_metrics_Occur$EtOH_Occur*100 == 0, 
                          df_enviro_metrics_Occur$EtOH_Occur*100 + 5, 
                          df_enviro_metrics_Occur$EtOH_Occur*100 - 5)), 
             colour = "steelblue", alpha=0.3, size = 2, position = position_jitter(height=5, width=0)) +
  ylim(0, 100) ; plot_pred_Occur_NDVI

#Finding model-predicted maximum values and confidence intervals
  #First name the variables (note we are writing over the value names from earlier):
x<-pred_Occur_NDVI$NDVI.mean #x (mean NDVI)
y<-pred_Occur_NDVI$pred_probOccur_NDVI #y (% occurrence)
upperCI<-pred_Occur_NDVI$pred_probOccur_NDVI_upper #upper confidence interval
lowerCI<-pred_Occur_NDVI$pred_probOccur_NDVI_lower #lower confidence interval
  #Calculate values:
max_index <- which.max(y) #index at maximum y
max_x <- x[max_index] ; max_x #value of x at max y (= 0.5)
max_y <- y[max_index]*100 ; max_y #max y (= 85.2% occurrence)
upperCI[max_index]*100 #upper CI (= 94.6%)
lowerCI[max_index]*100 #lower CI (= 65.6%)


####
#Code Segment 4.4) Environmental Analysis - Ethanol Concentration------------------

##Note: To run this analysis with a different minimum sampling threshold, return to Code Segment 4.2 Step 2, choose a different threshold and re-run Step 3 before returning here. While the results are reasonably consistent across sampling thresholds, it is worth noting that each threshold yields different "best" models depending on the combination of predictor variables considered, and we prioritized focusing on a single threshold and maximizing sample size. The code included below to calculate model predictions corresponds only to the most supported model based on this threshold, matching those reported in the manuscript.

#Like in Code Segment 4.3, We will build several Generalized Linear Mixed Models with our variables of interest, modeled as various combinations of linear and second-order polynomials and controlling for spatial autocorrelation, and then compare these models with corrected Akaike Information Criterion (AICc). We will do this both for mean and for maximum ethanol concentrations, and select the model that has the most support amongst both mean and maximum ethanol.

#First select records with positive (non-zero) ethanol measurements only
df_enviro_metrics_ConcPos <- df_enviro_metrics_Conc %>% filter(EtOH_mean>0)

#Visualize distribution of positive proportion mean ethanol values:
hist(df_enviro_metrics_ConcPos$EtOH_mean, xlab="Mean Ethanol (proportion)") 
#appears very left skewed and zero-inflated
#We will fit a Gamma distribution to each model, given the distribution shape


#Run models for mean ethanol
glmm_EtOHmean_1 <- glmmTMB(EtOH_mean ~ -1 + realm + scale(bio1) + scale(bio12) + scale(bio15) + scale(NDVI.mean) +  
                             exp(pos + 0|1), family = Gamma(link = "log"), data=df_enviro_metrics_ConcPos)

glmm_EtOHmean_2 <- glmmTMB(EtOH_mean ~ -1 + realm + poly(bio1,2) + poly(bio12,2) + poly(bio15,2) + poly(NDVI.mean,2) +  
                             exp(pos + 0|1), family = Gamma(link = "log"), data=df_enviro_metrics_ConcPos)

glmm_EtOHmean_3 <- glmmTMB(EtOH_mean ~ scale(bio1) + scale(bio12) + scale(bio15) + scale(NDVI.mean) +  
                             exp(pos + 0|1), family = Gamma(link = "log"), data=df_enviro_metrics_ConcPos)

glmm_EtOHmean_4 <- glmmTMB(EtOH_mean ~ poly(bio1,2) +  poly(bio12,2) +  poly(bio15,2) + poly(NDVI.mean,2) +  
                             exp(pos + 0|1), family = Gamma(link = "log"), data=df_enviro_metrics_ConcPos)

glmm_EtOHmean_5 <- glmmTMB(EtOH_mean ~ poly(bio1,2) +  poly(bio12,2) +  poly(bio15,2) +  
                             exp(pos + 0|1), family = Gamma(link = "log"), data=df_enviro_metrics_ConcPos)

glmm_EtOHmean_6 <- glmmTMB(EtOH_mean ~ -1 + realm +  
                             exp(pos + 0|1), family = Gamma(link = "log"), data=df_enviro_metrics_ConcPos)

glmm_EtOHmean_7 <- glmmTMB(EtOH_mean ~ bio1 +  
                             exp(pos + 0|1), family = Gamma(link = "log"), data=df_enviro_metrics_ConcPos)

glmm_EtOHmean_8 <- glmmTMB(EtOH_mean ~ -1 + realm + poly(bio1,2) + poly(bio12,2) +  
                             exp(pos + 0|1), family = Gamma(link = "log"), data=df_enviro_metrics_ConcPos)

glmm_EtOHmean_9 <- glmmTMB(EtOH_mean ~ -1 + realm + poly(bio1,2)  +  
                             exp(pos + 0|1), family = Gamma(link = "log"), data=df_enviro_metrics_ConcPos)

glmm_EtOHmean_10 <- glmmTMB(EtOH_mean ~ poly(bio1,2)  +  
                              exp(pos + 0|1) + (realm|1), family = Gamma(link = "log"), data=df_enviro_metrics_ConcPos)

#Extract AICc values for all models
models_AICc_mean <- MuMIn::AICc(glmm_EtOHmean_1, glmm_EtOHmean_2, glmm_EtOHmean_3, glmm_EtOHmean_4, glmm_EtOHmean_5, glmm_EtOHmean_6, glmm_EtOHmean_7, glmm_EtOHmean_8, glmm_EtOHmean_9, glmm_EtOHmean_10)
models_AICc_mean[order(models_AICc_mean$AICc), ] #order from lowest to highest


#Run models for maximum ethanol
glmm_EtOHmax_1 <- glmmTMB(EtOH_max ~ -1 + realm + scale(bio1) + scale(bio12) + scale(bio15) + scale(NDVI.mean) +  
                            exp(pos + 0|1), family = Gamma(link = "log"), data=df_enviro_metrics_ConcPos)

glmm_EtOHmax_2 <- glmmTMB(EtOH_max ~ -1 + realm + poly(bio1,2) + poly(bio12,2) + poly(bio15,2) + poly(NDVI.mean,2) +  
                            exp(pos + 0|1), family = Gamma(link = "log"), data=df_enviro_metrics_ConcPos)

glmm_EtOHmax_3 <- glmmTMB(EtOH_max ~ scale(bio1) + scale(bio12) + scale(bio15) + scale(NDVI.mean) +  
                            exp(pos + 0|1), family = Gamma(link = "log"), data=df_enviro_metrics_ConcPos)

glmm_EtOHmax_4 <- glmmTMB(EtOH_max ~ poly(bio1,2) +  poly(bio12,2) +  poly(bio15,2) + poly(NDVI.mean,2) +  
                            exp(pos + 0|1), family = Gamma(link = "log"), data=df_enviro_metrics_ConcPos)

glmm_EtOHmax_5 <- glmmTMB(EtOH_max ~ poly(bio1,2) +  poly(bio12,2) +  poly(bio15,2) +  
                            exp(pos + 0|1), family = Gamma(link = "log"), data=df_enviro_metrics_ConcPos)

glmm_EtOHmax_6 <- glmmTMB(EtOH_max ~ -1 + realm +  
                            exp(pos + 0|1), family = Gamma(link = "log"), data=df_enviro_metrics_ConcPos)

glmm_EtOHmax_7 <- glmmTMB(EtOH_max ~ bio1 +  
                            exp(pos + 0|1), family = Gamma(link = "log"), data=df_enviro_metrics_ConcPos)

glmm_EtOHmax_8 <- glmmTMB(EtOH_max ~ -1 + realm + poly(bio1,2) + poly(bio12,2) +  
                            exp(pos + 0|1), family = Gamma(link = "log"), data=df_enviro_metrics_ConcPos)

glmm_EtOHmax_9 <- glmmTMB(EtOH_max ~ -1 + realm + poly(bio1,2)  +  
                            exp(pos + 0|1), family = Gamma(link = "log"), data=df_enviro_metrics_ConcPos)

glmm_EtOHmax_10 <- glmmTMB(EtOH_max ~ poly(bio1,2)  +  
                             exp(pos + 0|1), family = Gamma(link = "log"), data=df_enviro_metrics_ConcPos)

#Extract AICc values for all models
models_AICc_max <- MuMIn::AICc(glmm_EtOHmax_1, glmm_EtOHmax_2, glmm_EtOHmax_3, glmm_EtOHmax_4, glmm_EtOHmax_5, glmm_EtOHmax_6, glmm_EtOHmax_7, glmm_EtOHmax_8, glmm_EtOHmax_9, glmm_EtOHmax_10)
models_AICc_max[order(models_AICc_max$AICc), ] #order from lowest to highest


#Model 5 is the most supported (i.e. lowest AICc) amongst *both* mean and maximum ethanol models
#Print results for both Model 5's
summary(glmm_EtOHmean_5)
summary(glmm_EtOHmax_5)
#BIO12 and BIO15 both have significant impact on ethanol concentration (p < 0.05)


##


#Calculate model predictions (for Model #5 for both mean and maximum ethanol)

#Assign alpha and z-scores so we can calculate confidence intervals on the log scale
alpha <- 0.05
z_score <- qnorm(1 - alpha / 2)  # 95% confidence interval


  #1) BIO12 (Annual Precipitation)

#First create dataframe of the predictor variables, filled with raw BIO12 values and mean values of all other predictors included in the model
pred_Conc_bio12 <- data.frame(
  bio1 = mean(df_enviro_metrics_ConcPos$bio1),
  bio12 = df_enviro_metrics_ConcPos$bio12,
  bio15 = mean(df_enviro_metrics_ConcPos$bio15),
  pos = df_enviro_metrics_ConcPos$pos)

#Calculate predictions for partial response BIO12
  #(a) Mean ethanol concentration
    #Predictions:
predictions_EtOHmean_bio12 <- predict(glmm_EtOHmean_5, pred_Conc_bio12, type = "response", se.fit = TRUE)
    #Calculate fit and upper/lower bounds of the predictions
pred_Conc_bio12$pred_EtOHmean <- predictions_EtOHmean_bio12$fit 
pred_Conc_bio12$pred_EtOHmean_lower <- predictions_EtOHmean_bio12$fit - z_score * predictions_EtOHmean_bio12$se.fit
pred_Conc_bio12$pred_EtOHmean_upper <- predictions_EtOHmean_bio12$fit + z_score * predictions_EtOHmean_bio12$se.fit
  #(b) Maximum ethanol concentration
    #Predictions:
predictions_EtOHmax_bio12 <- predict(glmm_EtOHmax_5, pred_Conc_bio12, type = "response", se.fit = TRUE)
    #Calculate fit and upper/lower bounds of the predictions
pred_Conc_bio12$pred_EtOHmax <- predictions_EtOHmax_bio12$fit 
pred_Conc_bio12$pred_EtOHmax_lower <- predictions_EtOHmax_bio12$fit - z_score * predictions_EtOHmax_bio12$se.fit
pred_Conc_bio12$pred_EtOHmax_upper <- predictions_EtOHmax_bio12$fit + z_score * predictions_EtOHmax_bio12$se.fit

#Plot the model predictions with 95% confidence intervals (Figure 4A)
plot_pred_Conc_bio12 <- ggplot(pred_Conc_bio12, aes(x = bio12)) +
  #Add mean EtOH predicted fit & confidence intervals:
  geom_line(aes(y = pred_EtOHmean*100, color = "Mean Ethanol"), lwd=2) +
  geom_ribbon(aes(ymin = pred_EtOHmean_lower*100, ymax = pred_EtOHmean_upper*100), fill = "grey80", alpha = 0.2) +
  #Add maximum EtOH predicted fit & confidence intervals:
  geom_line(aes(y = pred_EtOHmax*100, color = "Maximum Ethanol"), lwd=2) +
  geom_ribbon(aes(ymin = pred_EtOHmax_lower*100, ymax = pred_EtOHmax_upper*100), fill = "grey60", alpha = 0.2) +  
  #Further aesthetics:
  scale_color_manual(values = c("Mean Ethanol" = "deepskyblue1", "Maximum Ethanol" = "deepskyblue4")) + 
  labs(title = "", 
       x = "Annual Precipitation (mm)", 
       y = "Ethanol Concentration (%)",
       color = "Legend") +
  theme_bw() + 
  #Add a layer with raw data points for maximum ethanol:
  geom_point(aes(x=df_enviro_metrics_ConcPos$bio12, 
                 y=df_enviro_metrics_ConcPos$EtOH_max*100), 
             colour="deepskyblue4", alpha=0.4) +
  #Add a layer with raw data points for mean ethanol:
  geom_point(aes(x=df_enviro_metrics_ConcPos$bio12, 
                 y=df_enviro_metrics_ConcPos$EtOH_mean*100), 
             colour="deepskyblue1", alpha=0.3) ; plot_pred_Conc_bio12

#Finding model-predicted maximum values and confidence intervals
  #First name the variables:
x<-pred_Conc_bio12$bio12 #x (bio12)
ymean<-pred_Conc_bio12$pred_EtOHmean #y for mean EtOH concentration
ymax<-pred_Conc_bio12$pred_EtOHmax #y for maximum EtOH concentration
upperCImean<-pred_Conc_bio12$pred_EtOHmean_upper #upper confidence interval (for mean EtOH)
lowerCImean<-pred_Conc_bio12$pred_EtOHmean_lower #lower confidence interval (for mean EtOH)
upperCImax<-pred_Conc_bio12$pred_EtOHmax_upper #upper confidence interval (for maximum EtOH)
lowerCImax<-pred_Conc_bio12$pred_EtOHmax_lower #lower confidence interval (for maximum EtOH)
  #Calculate values:
    #(a) Mean ethanol concentration:
max_index_mean <- which.max(ymean) #index at maximum y
x[max_index_mean] #value of x at max y (= 3335 mm)
ymean[max_index_mean]*100 #max y value (= 0.46% ABV)
upperCImean[max_index_mean]*100 #upper CI (= 0.97%)
lowerCImean[max_index_mean]*100 #lower CI (= -0.05%, but since we cannot have negative % ABV, this effectively = 0)
    #(b) Maximum ethanol concentration:
max_index_max <- which.max(ymax) #index at maximum y
x[max_index_max] #value of x at max y (= 3335 mm)
ymax[max_index_max]*100 #max y value (= 1.39% ABV)
upperCImax[max_index_max]*100 #upper CI (= 2.83%)
lowerCImax[max_index_max]*100 #lower CI (= -0.04%, but since we cannot have negative % ABV, this effectively = 0)


  #2) BIO15 (Precipitation Seasonality)

#First create dataframe of the predictor variables, filled with raw BIO15 values and mean values of all other predictors included in the model
pred_Conc_bio15 <- data.frame(
  bio1 = mean(df_enviro_metrics_ConcPos$bio1),
  bio12 = mean(df_enviro_metrics_ConcPos$bio12),
  bio15 = df_enviro_metrics_ConcPos$bio15,
  pos = df_enviro_metrics_ConcPos$pos)

#Calculate predictions for partial response BIO12
  #(a) Mean ethanol concentration
    #Predictions:
predictions_EtOHmean_bio15 <- predict(glmm_EtOHmean_5, pred_Conc_bio15, type = "response", se.fit = TRUE)
    #Calculate fit and upper/lower bounds of the predictions
pred_Conc_bio15$pred_EtOHmean <- predictions_EtOHmean_bio15$fit 
pred_Conc_bio15$pred_EtOHmean_lower <- predictions_EtOHmean_bio15$fit - z_score * predictions_EtOHmean_bio15$se.fit
pred_Conc_bio15$pred_EtOHmean_upper <- predictions_EtOHmean_bio15$fit + z_score * predictions_EtOHmean_bio15$se.fit
  #(b) Maximum ethanol concentration
    #Predictions:
predictions_EtOHmax_bio15 <- predict(glmm_EtOHmax_5, pred_Conc_bio15, type = "response", se.fit = TRUE)
    #Calculate fit and upper/lower bounds of the predictions
pred_Conc_bio15$pred_EtOHmax <- predictions_EtOHmax_bio15$fit 
pred_Conc_bio15$pred_EtOHmax_lower <- predictions_EtOHmax_bio15$fit - z_score * predictions_EtOHmax_bio15$se.fit
pred_Conc_bio15$pred_EtOHmax_upper <- predictions_EtOHmax_bio15$fit + z_score * predictions_EtOHmax_bio15$se.fit

#Plot the model predictions with 95% confidence intervals (Figure 4B)
plot_pred_Conc_bio15 <- ggplot(pred_Conc_bio15, aes(x = bio15)) +
  #Add mean EtOH predicted fit & confidence intervals:
  geom_line(aes(y = pred_EtOHmean*100, color = "Mean Ethanol"), lwd=2) +
  geom_ribbon(aes(ymin = pred_EtOHmean_lower*100, ymax = pred_EtOHmean_upper*100), fill = "grey80", alpha = 0.2) +
  #Add maximum EtOH predicted fit & confidence intervals:
  geom_line(aes(y = pred_EtOHmax*100, color = "Maximum Ethanol"), lwd=2) +
  geom_ribbon(aes(ymin = pred_EtOHmax_lower*100, ymax = pred_EtOHmax_upper*100), fill = "grey60", alpha = 0.2) +  
  #Further aesthetics:
  scale_color_manual(values = c("Mean Ethanol" = "deepskyblue1", "Maximum Ethanol" = "deepskyblue4")) + 
  labs(title = "", 
       x = "Precipitation Seasonality (CV)", 
       y = "Ethanol Concentration (%)",
       color = "Legend") +
  theme_bw() + 
  #Add a layer with raw data points for maximum ethanol:
  geom_point(aes(x=df_enviro_metrics_ConcPos$bio15, 
                 y=df_enviro_metrics_ConcPos$EtOH_max*100), 
             colour="deepskyblue4", alpha=0.4) +
  #Add a layer with raw data points for mean ethanol:
  geom_point(aes(x=df_enviro_metrics_ConcPos$bio15, 
                 y=df_enviro_metrics_ConcPos$EtOH_mean*100), 
             colour="deepskyblue1", alpha=0.3) ; plot_pred_Conc_bio15

#Finding model-predicted maximum values and confidence intervals
  #First name the variables (note we are writing over the value names from earlier):
x<-pred_Conc_bio15$bio15 #x (bio15)
ymean<-pred_Conc_bio15$pred_EtOHmean #y for mean EtOH concentration
ymax<-pred_Conc_bio15$pred_EtOHmax #y for maximum EtOH concentration
upperCImean<-pred_Conc_bio15$pred_EtOHmean_upper #upper confidence interval (for mean EtOH)
lowerCImean<-pred_Conc_bio15$pred_EtOHmean_lower #lower confidence interval (for mean EtOH)
upperCImax<-pred_Conc_bio15$pred_EtOHmax_upper #upper confidence interval (for maximum EtOH)
lowerCImax<-pred_Conc_bio15$pred_EtOHmax_lower #lower confidence interval (for maximum EtOH)
  #Calculate values:
    #(a) Mean ethanol concentration:
max_index_mean <- which.max(ymean) #index at maximum y
x[max_index_mean] #value of x at max y (= 91.7% mm)
ymean[max_index_mean]*100 #max y value (= 0.17% ABV)
upperCImean[max_index_mean]*100 #upper CI (= 0.23%)
lowerCImean[max_index_mean]*100 #lower CI (= 0.12%)
    #(b) Maximum ethanol concentration:
max_index_max <- which.max(ymax) #index at maximum y
x[max_index_max] #value of x at max y (= 91.7%)
ymax[max_index_max]*100 #max y value (= 0.54% ABV)
upperCImax[max_index_max]*100 #upper CI (= 0.71%)
lowerCImax[max_index_max]*100 #lower CI (= 0.36%)


####
#Code Segment 5) Canopy vs Ground Analysis--------------

  #Step 1) Filter the repeatedly sampled dataset further to obtain species which have at least 1 sample from the canopy ("plant") and from the forest floor ("ground"), AND had ethanol detected in at least one sample across entire raw data

#First select the sampling threshold to use going forward (these are obtained from Code Segment 1); for final analysis, we used Option 3 
  #Option 1) species with at least 3 ripe fruits: 
    #run code:     
      #df_repeated_samples <- df_3R
  #Option 2) species with at least 3 ripe or overripe fruits:
    #run code:     
      #df_repeated_samples <- df_3ROR
  #Option 3) species with at least 3 midripe, ripe, or overripe fruits: 
    df_repeated_samples <- df_3MRROR

#Apply further filtering:
#Note we need to remove records retained in the full fruit dataset which do not have a reliable % ethanol value (i.e. records maintained solely as evidence of ethanol presence/absence)
df_GvsP <- df_repeated_samples %>%
  filter(ExcludeFromFinalAnalysis=="No") %>% #include only records with reliable %-EtOH values
  filter(EthanolDetectedInSpecies == "Yes") %>% #include only species with ethanol detected
  filter(LocationDescription=="ground" | LocationDescription=="plant") %>% #include only 'ground' and 'plant' collection sites (thus excluding 'unknown' or other values)
  group_by(FruitSpeciesName) %>%
  filter(any(LocationDescription=="ground") & any(LocationDescription=="plant")) %>% #include only species that have at least 1 plant and 1 ground sample
  ungroup()

#Assign relevant variables as factors
df_GvsP$LocationDescription <- as.factor(df_GvsP$LocationDescription)
levels(df_GvsP$LocationDescription) #check levels -- ground is first level
df_GvsP$IntactWhenMeasured <- as.factor(df_GvsP$IntactWhenMeasured)


  #Step 2) Investigate distribution of the dependent variable (Proportion EtOH), including normality and variance of residuals.

#Visualize distribution of proportion ethanol:
hist(df_GvsP$ProportionEtOH, xlab="Ethanol Produced (proportion)") 
#appears very left skewed and zero-inflated

#Normality:
  #Visual test with QQ-plot
qqnorm(df_GvsP$ProportionEtOH) 
#does not appear normally distributed

  #Formal test with Shapiro-Wilk normality test
shapiro.test(df_GvsP$ProportionEtOH)
#not normally distributed (p < 0.05)

#Compare sample sizes between groups 
length(df_GvsP$ProportionEtOH[df_GvsP$LocationDescription=="ground"]) #1240 measurements
length(df_GvsP$ProportionEtOH[df_GvsP$LocationDescription=="plant"]) #811 measurements
#unbalanced sample sizes

#Formal test of homoscedasticity:
  #Use Fligner-Killeen test, which is robust against normality violations
fligner.test(df_GvsP$ProportionEtOH ~ df_GvsP$LocationDescription)
#residuals are not homoscedastic (p < 0.05)


  #Step 3) Modelling

#We will build Generalized Linear Mixed Models with variables of interest, both with and without zero-inflation, and compare models with corrected Akaike Information Criteria (AICc) values. We will fit a beta distribution given the data type (proportion) and distribution shape of the dependent variable.

#Build the "conditional" model (i.e. no zero-inflation)

  #First transform y-variable (ProportionEtOH) using Smithson and Verkuilen transformation so there are no longer exact 0s and we can fit a beta distribution (which is constricted to values between 0 and 1)
y <- df_GvsP$ProportionEtOH
n <- sum(!is.na(y))
ProportionEtOH_trans <- (y * (n - 1) + 0.5) / n

  #Run model
GvsP_conditional <- glmmTMB(ProportionEtOH_trans ~ LocationDescription + (1 | FruitSpeciesName) + (1 | IntactWhenMeasured) + (1|FruitUniqueID) + (1|SampleUniqueID),
                                      data = df_GvsP,
                                      ziformula = ~0, #0 means no zero-inflation
                                      family=beta_family())
summary(GvsP_conditional) #print model results


#Next build the same model controlling for zero-inflation
GvsP_zi <- glmmTMB(ProportionEtOH ~ LocationDescription + (1 | FruitSpeciesName) + (1 | IntactWhenMeasured) + (1|FruitUniqueID) + (1|SampleUniqueID),
                            data = df_GvsP,
                            ziformula = ~., #~. means the same as the conditional formula
                            family=beta_family())
summary(GvsP_zi) #print model results
#results are qualitatively fairly similar between models


#Compare models with AICc:
GvsP_models_AICc <- MuMIn::AICc(GvsP_conditional, GvsP_zi)
GvsP_models_AICc[order(GvsP_models_AICc$AICc), ] #order AICcs from lowest to highest
#conditional (i.e. non-zero-inflated model) is the most supported model (i.e. lowest AICc)


  #Step 4) Extract estimated marginal means (EMMs) from the best model

GvsP_EMMs <- emmeans(GvsP_conditional, ~LocationDescription, component="response") #extract EMMs
GvsP_EMMs_df <- as.data.frame(GvsP_EMMs) #convert to dataframe

#Transform EMMs from proportion to percentage:
GvsP_EMMs_percentage <- data.frame(
  LocationDescription = GvsP_EMMs_df$LocationDescription,
  emmean_percentage = GvsP_EMMs_df$emmean*100,
  SE_percentage = GvsP_EMMs_df$SE*100,
  df = GvsP_EMMs_df$df,
  LCL_percentage = GvsP_EMMs_df$asymp.LCL*100,
  UCL_percentage = GvsP_EMMs_df$asymp.UCL*100)
print(GvsP_EMMs_percentage)

#Calculate one-sided pairwise comparison (z-test) of the EMMs
  #Note our alternative hypothesis is that ethanol in ground-collected fruits is greater than ethanol in plant-collected fruits; our null hypothesis is that there is no difference in ethanol between ground- and plant-collected fruits:
GvsP_EMMs_pairwise <- pairs(GvsP_EMMs, side=">") 
GvsP_EMMs_pairwise #print results
#p < 0.05 --> reject null hypothesis and accept alternative hypothesis

#Transform pairwise EMMs from proportion to percentage (as is reported in manuscript)
GvsP_EMMs_pairwise_df <- as.data.frame(GvsP_EMMs_pairwise) #write as data frame
GvsP_EMMs_pairwise_percentage <- data.frame(
  Estimate.percentage = GvsP_EMMs_pairwise_df$estimate*100,
  SE.percentage = GvsP_EMMs_pairwise_df$SE*100,
  df = GvsP_EMMs_pairwise_df$df,
  z.ratio.percentage = GvsP_EMMs_pairwise_df$z.ratio*100,
  p.value = GvsP_EMMs_pairwise_df$p.value)
print(GvsP_EMMs_pairwise_percentage)

#Plot ethanol distribution for fruits collected from canopy versus ground (Figure 5)
location_order <- factor(df_GvsP$LocationDescription, level = c("plant", "ground")) #re-order ground and canopy for the plot
ripeness_order <- factor(df_GvsP$RipenessStage, level = c("midripe", "ripe", "overripe")) #re-order ripeness stage for the plot
plot_GvsP <- ggplot(data = df_GvsP, aes(x = location_order, y = ProportionEtOH*100, color = ripeness_order)) +
  #Add raw data points:
  geom_point(alpha=0.4, 
             size = 3,
             position = position_jitter(0.1)) +
  #Add layer containing error bars:
  geom_errorbar(data = GvsP_EMMs_percentage, 
                aes(x= LocationDescription, y = emmean_percentage,
                    ymin = emmean_percentage - SE_percentage, 
                    ymax = emmean_percentage + SE_percentage), 
                linewidth = 1.5, width = 0.075, colour = "black", inherit.aes = FALSE) +
  #Add layer containing EMMs:
  geom_point(data = GvsP_EMMs_percentage,
             aes(x = LocationDescription, y = emmean_percentage),
             size = 2, colour = "black", inherit.aes = FALSE) +
  #Further aesthetics:
  scale_color_manual(values = c("limegreen", "orange", "violetred4")) +
  labs(colour="Ripeness Stage") +
  ylab("Ethanol Concentration (%)") +  xlab("Sampling Location") +
  scale_x_discrete(labels=c('Canopy','Ground')) +
  theme_bw() ; plot_GvsP


  #Step 5) Run model with fruit ripeness stage as an interacting predictor variable

##Note we only ran this model including ripeness stage as a predictor variable on the least restrictive sampling threshold (i.e. at least three midripe, ripe, overripe fruits sampled per species), because this was the only threshold that included all ripeness stages (other than unripe).

#Check multicollinearity of RipenessStage and LocationDescription:

  #Test if the distribution of ripeness ranks differs significantly between canopy vs ground using a one-sided Wilcoxon rank sum test (hypothesis is that ground fruits are more ripe than canopy fruits)
wilcox.test(
  as.numeric(factor(RipenessStage, levels = c("midripe", "ripe", "overripe"), ordered = TRUE)) ~ LocationDescription, 
  alternative = "greater", 
  data = df_GvsP)
#p < 0.05, meaning ground-collected fruits are significantly more ripe than plant-collected fruits

  #Measure correlation strength between RipenessStage and LocationDescription using rank-biserial method
rank_biserial(
  as.numeric(factor(RipenessStage, levels = c("midripe", "ripe", "overripe"), ordered = TRUE)) ~ LocationDescription, 
  data = df_GvsP)
#r = 0.14, consistent with a weak positive correlation between the two variables. Given that the correlation is weak (i.e. canopy location is not a major driving predictor of how ripe the fruit is), we will continue forward with adding ripeness stage as a predictor variable in the model 


#Build the models to compare:
df_GvsP$RipenessStage <- factor(df_GvsP$RipenessStage, levels = c("midripe", "ripe", "overripe")) #re-level ripeness stage from least to most ripe

  #Conditional (i.e. non-zero-inflated) model:
#We will use the Smithson and Verkuilen transformed y-variable ('ProportionEtOH_trans') to enable fitting a beta distribution
GvsP_interaction_conditional <- glmmTMB(ProportionEtOH_trans ~ LocationDescription*RipenessStage + (1|FruitSpeciesName) + (1|IntactWhenMeasured) + (1|FruitUniqueID) + (1|SampleUniqueID),
                                data = df_GvsP,
                                ziformula = ~0,
                                family=beta_family())
summary(GvsP_interaction_conditional) #print model results

  #Zero-inflated model:
GvsP_interaction_zi <- glmmTMB(ProportionEtOH ~ LocationDescription*RipenessStage + (1|FruitSpeciesName) + (1|IntactWhenMeasured) + (1|FruitUniqueID) + (1|SampleUniqueID),
                                        data = df_GvsP,
                                        ziformula = ~.,
                                        family=beta_family())
summary(GvsP_interaction_zi) #print model results
#results are qualitatively fairly similar between models

  #Compare models with corrected Aikaike information criterion (AICc):
GvsP_interaction_models_AICc <- MuMIn::AICc(GvsP_interaction_conditional, GvsP_interaction_zi)
GvsP_interaction_models_AICc[order(GvsP_interaction_models_AICc$AICc), ] #order AICcs from lowest to highest
#conditional (i.e. non-zero-inflated model) has is the most supported model (i.e. lowest AICc)


#Run Type III Wald chi-square test on the best model
options(contrasts=c("contr.sum", "contr.poly")) #command ensures proper defaults for encoding variables
Anova(GvsP_interaction_conditional, type="III")
#the interaction term is significant (p < 0.05), so we will do some post-hoc contrasts


  #Step 6) Post-hoc comparisons to flesh out the significant interaction

#Extract estimated marginal means from the most supported model:
interaction_EMMs <- emmeans(GvsP_interaction_conditional, ~LocationDescription*RipenessStage)
interaction_EMMs #check the order of grouped variables to set up contrasts

#Write the contrasts accordingly:
interactioncontrasts <- list(
  C1=c(1,-1,0,0,0,0), #comparing ground vs plant for midripe fruits
  C2=c(0,0,1,-1,0,0), #comparing ground vs plant for ripe fruits
  C3=c(0,0,0,0,1,-1), #comparing ground vs plant for overripe fruits
  C4=c(1,0,-1,0,0,0), #comparing ground midripe vs ground ripe
  C5=c(0,0,1,0,-1,0), #comparing ground ripe vs ground overripe
  C6=c(0,1,0,-1,0,0), #comparing plant midripe vs plant ripe
  C7=c(0,0,0,1,0,-1) #comparing plant ripe vs plant overripe
)

#Run all pairwise comparisons using a Holm-Bonferonni p-value adjustment, and back-transform values to the response scale
#Note we will use two-sided, not one-sided, contrasts since these tests are post-hoc 
contrast(interaction_EMMs, method = interactioncontrasts, adjust="holm", type="response")
#C1 --> p > 0.05, NOT significant (G=P for midripe fruits)
#C2 --> p < 0.05, significant (G>P for ripe fruits)
#C3 --> p > 0.05, NOT significant (G=P for overripe fruits)
#C4 --> p < 0.05, significant (R>MR for ground fruits)
#C5 --> p < 0.05, significant (OR>R for ground fruits)
#C6 --> p > 0.05, NOT significant (MR=R for plant fruits)
#C7 --> p < 0.05, significant (OR>R for plant fruits)

#Because C7 is significant but C6 is not, we'll re-evaluate the contrast set now with an additional contrast:
interactioncontrasts_eight <- list(
  C1=c(1,-1,0,0,0,0), 
  C2=c(0,0,1,-1,0,0), 
  C3=c(0,0,0,0,1,-1), 
  C4=c(1,0,-1,0,0,0), 
  C5=c(0,0,1,0,-1,0), 
  C6=c(0,1,0,-1,0,0), 
  C7=c(0,0,0,1,0,-1), 
  C8=c(0,1/2,0,1/2,0,-1) #comparing average of plant ripe & plant midripe, vs plant overripe
)
contrast(interaction_EMMs, method = interactioncontrasts_eight, adjust="holm", type="response")
#no change in result interpretations for C1 thru C7
#C8 --> p < 0.05, significant (OR>MR=R for plant fruits)

#Extract 95% confidence intervals for all comparisons (as reported in manuscript)
confint(contrast(interaction_EMMs, method = interactioncontrasts_eight, adjust="holm", type="response"))


###
#Code Segment 6) Estimating Upper Ranges of Dietary Ethanol Exposure in Primates--------------

  #Step 1) Calculate plant species-specific ethanol metrics
  
##Note we use the Full Fruit Dataset (as opposed to the repeatedly sampled dataset), excluding unripe fruits, because this maximized the number of plant species with overlapping available ethanol measurements & primate foraging data. Since we are only focusing on upper ranges of ethanol exposure, the issue of using repeated sampling to reduce the chance of false negatives (i.e. species for which we did not detect ethanol but may have found ethanol if we sampled more fruits) is not as relevant, and we chose to prioritize sample size.

##Note the Full Fruit Dataset ('df_full_fruit') was obtained in Code Segment 1

##Also note we will calculate plant species-specific metrics that are specific to the sampling locations (i.e. field sites) for which we also have primate foraging data.

df_fruit_exposure <- df_full_fruit %>%
  filter(ExcludeFromFinalAnalysis=="No") %>% #filter for only records for which we have reliable %-EtOH values
  filter(RipenessStage=="midripe" | RipenessStage=="ripe" | RipenessStage=="overripe") %>% #filter for only ripened fruits (i.e. no unripe)
  filter(StudySiteName=="Phu Khieo Wildlife Sanctuary" | StudySiteName=="Kianjavato" | StudySiteName=="Tuanan" | StudySiteName=="Kibale" | StudySiteName=="Sector Santa Rosa" | StudySiteName=="Yakushima Island") %>% #filter for only the 6 field sites for which we have primate foraging data
  mutate(species_site_ID = paste0(FruitSpeciesName, "_", StudySiteName)) %>% #create a new column that combines fruit species name and field site name
  group_by(species_site_ID) %>%
  summarise(mean_EtOH_prop = mean(ProportionEtOH), #mean EtOH (proportion)
            max_EtOH_prop = max(ProportionEtOH)) #maximum EtOH (proportion)


#Combine these metrics with the primate foraging data

df_primate_exposure <- read.csv("Casorso_et_al_2026_-_Primate_Foraging_Data.csv") #load primate foraging data

df_primate_exposure <- df_primate_exposure %>%
  mutate(species_site_ID = paste0(Plant_species, "_", Site_Name)) #create a new column that combines fruit species name and field site name

df_exposure <- left_join(df_primate_exposure, df_fruit_exposure, by="species_site_ID") #combine the two dataframes together
df_exposure <- df_exposure %>% filter(!is.na(mean_EtOH_prop)) #remove foraging records for plant species that do not have site-specific ethanol measurements

#Assign new column names for ease of coding
df_exposure <- df_exposure %>%
  rename(IntakeRate = Intake_rate_.g_wet_weight_per_min.,
         FruitMassConsumed = Fruit_consumed_.g.wet.weight.,
         ForagingDuration = Foraging_duration_.min.,
         BodyWeight = Average_Adult_Body_Weight_.kg.,
         Frugivory = Average_Percent_Frugivory)

#Ensure variables are formatted properly
df_exposure$FruitMassConsumed <- as.numeric(df_exposure$FruitMassConsumed, na.rm=TRUE)
df_exposure$IntakeRate <- as.numeric(df_exposure$IntakeRate, na.rm=TRUE)
df_exposure$ForagingDuration <- as.numeric(df_exposure$ForagingDuration, na.rm=TRUE)
df_exposure$Unique_Daily_ID <- as.factor(df_exposure$Unique_Daily_ID)
df_exposure$FruitMassConsumed <- as.numeric(df_exposure$FruitMassConsumed, na.rm=TRUE)
df_exposure$BodyWeight <- as.numeric(df_exposure$BodyWeight, na.rm=TRUE)
df_exposure$Frugivory <- as.numeric(df_exposure$Frugivory, na.rm=TRUE)


  #Step 2) Estimate ethanol intake

#(a) Calculate fruit wet weight consumed (grams wet weight) by multiplying foraging duration (min) by intake rate (grams wet weight / min) (intake rates are plant- and primate-species-specific)
df_exposure <- df_exposure %>%
  mutate(FruitMassConsumed = coalesce(FruitMassConsumed, ForagingDuration*IntakeRate)) #using the coalesce() function to calculate this for blank cells only (Kibale dataset already has mass wet weight consumed values calculated)


#(b) Calculate mass ethanol consumed (grams EtOH) by multiplying fruit mass consumed (grams wet weight) by the mean or maximum ethanol concentration (proportion EtOH)
df_exposure <- df_exposure %>%
  mutate(EtOH_consumed_mean = FruitMassConsumed*mean_EtOH_prop,
         EtOH_Consumed_max = FruitMassConsumed*max_EtOH_prop)


#(c) Obtain daily estimates of mass ethanol consumed (grams EtOH / day) by adding together any multiple values of mass ethanol consumed (grams EtOH) across a day for the same individual
  ##Note only the top 10 summed values of daily ethanol consumed (based on mean EtOH) and the top 10 summed values of daily ethanol consumed (based on maximum EtOH) were selected and presented in the manuscript. The supplied primate foraging data only contains the foraging records that yielded these top 10 mean-based values and top 10 max-based values. 

#We'll sum and select the top 10 individual daily consumption estimates separately for the mean-EtOH-based and max-EtOH-based values.

  #Values based on mean ethanol:
df_exposure_meanEtOH <- df_exposure %>%
  group_by(Primate_species, Unique_Daily_ID, BodyWeight, Frugivory) %>%
  summarise(Daily_EtOH_Consumed_mean = sum(EtOH_consumed_mean),
            .groups="drop") %>%
  group_by(Primate_species) %>% slice_max(Daily_EtOH_Consumed_mean, n=10) %>%
  ungroup()

  #Values based on  maximum ethanol:
df_exposure_maxEtOH <- df_exposure %>%
  group_by(Primate_species, Unique_Daily_ID, BodyWeight, Frugivory) %>%
  summarise(Daily_EtOH_Consumed_max = sum(EtOH_Consumed_max),
            .groups="drop") %>%
  group_by(Primate_species) %>% slice_max(Daily_EtOH_Consumed_max, n=10) %>%
  ungroup()

##Note that not every primate species had 10 mean-based or max-based values to choose from, due to a lack of available data (constrained by requiring both site-specific ethanol measurements and foraging data for the same plant species); see dataframes for which primate species this applied to.


#(d) Calculate daily ethanol dosages standardized by body mass (grams EtOH / kilogram / day), by dividing each daily ethanol consumption (grams EtOH / day) by the average adult body weight (kilograms) for that primate species

  #Values based on mean ethanol:
df_exposure_meanEtOH <- df_exposure_meanEtOH %>%
  mutate(DailyDosage = Daily_EtOH_Consumed_mean / BodyWeight)

  #Values based on maximum ethanol:
df_exposure_maxEtOH <- df_exposure_maxEtOH %>%
  mutate(DailyDosage = Daily_EtOH_Consumed_max / BodyWeight)


#(e) Calculate number daily standard drinks (# / day), by dividing daily dosages by a "human conversion factor" of 10/70 (i.e. 10-g EtOH in a standard drink for a 70-kg human); see manuscript for more details.

  #Values based on mean ethanol:
df_exposure_meanEtOH <- df_exposure_meanEtOH %>%
  mutate(DailyStandardDrinks = DailyDosage / (10/70))

  #Values based on maximum ethanol:
df_exposure_maxEtOH <- df_exposure_maxEtOH %>%
  mutate(DailyStandardDrinks = DailyDosage / (10/70))


#(f) For each primate species, calculate the mean and standard deviation for Daily Ehanol Consumption (g/day), Daily Ethanol Dosage (g/kg/day), and Daily Number Standard Drinks (#/day) (as is presented in the manuscript)

  #Values based on mean ethanol (Supplementary Table S8 (first portion)):
df_meanexposure_meanEtOH <- df_exposure_meanEtOH %>%
  group_by(Primate_species) %>%
  reframe(Mean_DailyEtOHConsumed = mean(Daily_EtOH_Consumed_mean),
          SD_DailyEtOHConsumed = sd(Daily_EtOH_Consumed_mean),
          Avg_DailyDosage = mean(DailyDosage),
          SD_DailyDosage = sd(DailyDosage),
          Avg_StandardDrinks = mean(DailyStandardDrinks),
          SD_StandardDrinks = sd(DailyStandardDrinks))
print(df_meanexposure_meanEtOH) #print results table

  #Values based on maximum ethanol (Supplementary Table S8 (second portion)):
df_meanexposure_maxEtOH <- df_exposure_maxEtOH %>%
  group_by(Primate_species) %>%
  reframe(Mean_DailyEtOHConsumed = mean(Daily_EtOH_Consumed_max),
          SD_DailyEtOHConsumed = sd(Daily_EtOH_Consumed_max),
          Avg_DailyDosage = mean(DailyDosage),
          SD_DailyDosage = sd(DailyDosage),
          Avg_StandardDrinks = mean(DailyStandardDrinks),
          SD_StandardDrinks = sd(DailyStandardDrinks))
print(df_meanexposure_maxEtOH) #print results table


#(g) Plotting upper range estimates of dietary ethanol exposure across primate species

  #Plot based on mean ethanol (Figure 6):
plot_mean <- ggplot(data = df_exposure_meanEtOH, aes(x = reorder(Primate_species, DailyStandardDrinks, FUN=mean), y = DailyStandardDrinks, fill=Frugivory, alpha=0.1, color=Frugivory)) + 
  geom_boxplot(outlier.shape=NA, 
               linewidth=0.8, median.linewidth = 0.7) + geom_point(pch=21) +  
  #Add a mean point to each boxplot:
  stat_summary(fun = mean, geom = "point", shape = 4, size = 4, alpha=0.9) +
  #Further aesthetics:
  xlab("Primate Species") + ylab("Standard Drinks (#/day)") + 
  scale_fill_gradient(high="purple", low="orange") + scale_color_gradient(high="purple",low="orange") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size=9, face='italic')) +
  labs(colour = "Frugivory (%)") +
  guides(fill="none", alpha="none") ; plot_mean


  #Plot based on maximum ethanol (Supplementary Figure S5):
plot_max <- ggplot(data = df_exposure_maxEtOH, aes(x = reorder(Primate_species, DailyStandardDrinks, FUN=mean), y = DailyStandardDrinks, fill=Frugivory, alpha=0.1, color=Frugivory)) + 
  geom_boxplot(outlier.shape=NA, 
               linewidth=0.8, median.linewidth = 0.7) + geom_point(pch=21) +  
  #Add a mean point to each boxplot:
  stat_summary(fun = mean, geom = "point", shape = 4, size = 4, alpha=0.9) +
  #Further aesthetics:
  xlab("Primate Species") + ylab("Standard Drinks (#/day)") + 
  scale_fill_gradient(high="purple", low="orange") + scale_color_gradient(high="purple",low="orange") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size=9, face='italic')) +
  labs(colour = "Frugivory (%)") +
  guides(fill="none", alpha="none") ; plot_max


###