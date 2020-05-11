# Yeast_Ionome_Metabolome

This repository contains materials for the paper [**Extraction and integration of genetic networks from short-profile omic datasets**], by J. Iacovacci et al (under review).

## Scripts
* *cosM.m* : function to measure the Mahalanobis cosine between the N observations of M correlated variables in a dataset.

* *cosH.m* : function to measure the hybrid-Mahalanobis cosine between the N observations of M correlated variables in a dataset.

* *multilink_analysis.m* : script to compute a basic multilink analysis of the Yeast Ionome-Metabolome multiplex network.


## Networks
The following files included in the repository contain the Yeast Ionome-Metabolome multiplex network:

* *yeast_ionome_metabolome_multiplex_cosM.txt* : edge lists of the layers of the Yeast Ionome-Metabolome multiplex network (ionome knockout, ionome oevrexpression, metabolome amino acids) extracted by using the Mahalanobis cosine. Link weights correspond to the rank of the correlation values. 
  
* *yeast_ionome_metabolome_multiplex_cosH.txt* : edge lists of the layers of the Yeast Ionome-Metabolome multiplex network (ionome knockout, ionome oevrexpression, metabolome amino acids) extracted by using the hybrid-Mahalanobis cosine. Link weights correspond to the rank of the correlation values. 
  
* *yeast_ionome_metabolome_multiplex_cos.txt* : edge lists of the layers of the Yeast Ionome-Metabolome multiplex network (ionome knockout, ionome oevrexpression, metabolome amino acids) extracted by using the cosine. Link weights correspond to the rank of the correlation values. 

* *yeast_ionome_metabolome_multiplex_PCC.txt* : edge lists of the layers of the Yeast Ionome-Metabolome multiplex network (ionome knockout, ionome oevrexpression, metabolome amino acids) extracted by using the Pearson's coefficient. Link weights correspond to the rank of the correlation values. 
    
* *nodes_list.txt* : contains the gene names (ORFs) corresponding to the nodes in the multiplex networks. 
  
 
 ## Data
The following data folders and files included in the repository contain datasets:

* *Synthetic_Dataset* : contains the files correspondng to the synthetic datasets generated for the analysis. 
  
* *profiles_all_datasets.mat* : contains all datasets in matlab format used to genearete the multiplex network (to be released upon manuscript acceptance).
  
  
  
