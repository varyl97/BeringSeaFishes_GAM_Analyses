###Generalized Additive Analyses: Pacific cod
#the following code creates generalized additive models for larvae of Pacific cod. 
#Pacific cod spawn beneath sea ice each year, therefore their eggs are rarely caught by ecoFOCI trawls. 
#   thus, Pac cod are only represented in the larval biogeography portion of my analyses. 

# Load larval data --------------------------------------------------------


apsub<-read.csv(file='../Ichthyo Data/Cleaned_Cut_ApEggs.csv',header=TRUE,check.names=TRUE)

aplarv.ctd<-read.csv(file='../Ichthyo Data/Cleaned_Cut_ApLarv_wCTD.csv',header=TRUE,check.names=TRUE)