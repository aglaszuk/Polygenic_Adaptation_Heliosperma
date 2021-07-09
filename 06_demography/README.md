
# Testing alternative demographic scenarios of ecotype pair divergence

We evaluated which demographic scenario better explains our data using the continuous-time coalescent simulator fastSimcoal2 v.2.6.0.3 (Excoffier et al. 2013). 
For each combination of two ecotype pairs we evaluated six models:

* SI: 1-origin
* SI: 2-origins
* IM: 1-origin, T2>T3
* IM: 1-origin, T3<T2
* IM: 2-origins, T2>T3
* IM: 2-origins, T3>T2

Where SI = "strict isolation", IM = "isolation with migration", and T2 and T3 represent the second and third split events in the model (see also Fig. 3 of the preprint).

Three files were required as input: the joint site frequency spectra (sfs) for all six combinations of four populations, a template file defining the demographic model and a estimation file defining the parameters priors and model rules. 
  
Examples input files for this analysis can be found in the input_files folder.
  
To launch several replicate runs of the anaysis on the Life Science Compute Cluster (LISC) (https://cube.univie.ac.at/lisc) we used the 00_launchFSC_HELIO.sh script.
The script 01_extract_bestL.sh was used to extract best likelihood average estimates across all runs and obtain the top 10 models parameters estimates to produce simulated sfs to be compared with the observed sfs.
