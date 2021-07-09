
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
