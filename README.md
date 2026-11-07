# Ethanol_in_Fruits_and_Flowers

## Data and code supporting Casorso, Wells, et al., titled _"Pervasive ethanol in the fruits and floral nectars consumed by tropical frugivores"_

Our repository contains one R script, 6 corresponding data files, and 1 corresponding metadata file. 


The R script **_“Casorso_et_al_2026_-_R Code.R”_** allows reproduction of all analyses included in the manuscript. It also allows reproduction of analytical comparison between each of the three minimum sampling thresholds we investigated for the Phylogenetic, Environmental, and Canopy vs. Ground analyses. Furthermore, code to reproduce the main figures and tables is included (aside from the phylogenetic trees, which were visualized using [iTOL](https://itol.embl.de/)).
The script contains the following segments, intended to be run in order:
1.  <ins>Code Segment 0)</ins> Preparation
    -	Installing/loading required packages
2.  <ins>Code Segment 1)</ins> Dataset Filtering
    - 	Requires file **_“Casorso_et_al_2026_-_raw_dataset.csv”_**
    -	Applies filters to raw dataset to obtain “full fruit dataset” and “repeatedly sampled dataset” used for subsequent analyses
3. <ins>Code Segment 2)</ins> Summary Statistics of Fruits & Floral Tissues
    - 	Calculates various descriptive statistics reported in the manuscript for fruits and floral tissue samples
3.  <ins>Code Segment 3.1)</ins> Phylogenetic Tree Creation
    - 	Requires the following files:
        -	 **_"Casorso_et_al_2026_-_Phylo_Tree_Species_List_3R.csv"_**
        -	 **_"Casorso_et_al_2026_-_Phylo_Tree_Species_List_3ROR.csv"_**
        -	 **_"Casorso_et_al_2026_-_Phylo_Tree_Species_List_3MRROR.csv"_**
        -	 **_"Casorso_et_al_2026_-_Phylo_Tree_Species_List_All.csv"_**
    - 	Creates phylogenetic tree of fruit species
4.  <ins>Code Segment 3.2)</ins> Phylogenetic Analysis (Pagel’s Lambda)
    -	 Calculates phylogenetic signal in median, mean, and maximum ethanol concentration using the tree obtained from Code Segment 3.1
5.  <ins>Code Segment 4.1)</ins> Extracting Environmental Data
    - 	Obtain publicly available environmental data to be used for analyses
    -   Requires map of realms/regions to be downloaded externally
        -  a) Follow link <https://macroecology.ku.dk/resources/wallace>.
        -	 b) Select the 'Download all maps' tab under 'GIS maps of CMEC zoogeographic realms and regions'.
        -	 c) Once the .zip file has downloaded, extract all the files to your chosen working directory to then be loaded into R.
6.  <ins>Code Segment 4.2)</ins> Binding Environmental Data to the Ethanol Measurements
    - 	Binds extracted environmental data from Code Segment 4.1 to the ethanol measurements
7.	<ins>Code Segment 4.3)</ins> Environmental Analysis – Ethanol Occurrence
    - 	Run environmental analysis to evaluate influence of predictor variables on ethanol occurrence, using bound data from Code Segment 4.2
8.	<ins>Code Segment 4.4)</ins> Environmental Analysis – Ethanol Concentration
    - 	Run environmental analysis to evaluate influence of predictor variables on ethanol concentration (mean and maximum), using bound data from Code Segment 4.2
9.	<ins>Code Segment 5)</ins> Canopy vs. Ground Analysis
    - 	Run canopy vs. ground analysis to evaluate influence of collection site and ripeness stage on ethanol concentration
10. <ins>Code Segment 6)</ins> Estimating Upper Ranges of Dietary Ethanol Exposure in Primates
     - 	Requires the file **_“Casorso_et_al_2026_-_Primate_Foraging_Data.csv”_**
     - 	Combines primate foraging data with ethanol measurements to estimate maximal daily consumption

If you want to run one analysis section independently out of order, you will need to first run Code Segments 0 and 1. Note Sections 3.1 to 3.2 and Sections 4.1 to 4.4 each <ins>_must_</ins> be run together in order.
Code is intended to be run in R version 4.6.0 (or later).



The metadata file **_“Casorso_et_al_2026_-_Metadata.xlsx”_** contains six worksheets, each of which corresponds to one of the six datafiles. Each worksheet contains information about the purpose and contents of the associated datafile, as well as a table with a description of each column heading.
