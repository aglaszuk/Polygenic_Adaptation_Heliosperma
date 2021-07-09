
# Testing alternative demographic scenarios of ecotype pair divergence

We evaluated which demographic scenario better explains our data using the continuous-time coalescent simulator fastSimcoal2 v.2.6.0.3 (Excoffier et al. 2013). 
For each combination of two ecotype pairs we evaluated six models:

* SI: 1-origin
* SI: 2-origins
* IM: 1-origin, T2>T3
* IM: 1-origin, T3<T2
* IM: 2-origins, T2>T3
* IM: 2-origins, T3>T2

Where SI = "strict isolation" and IM = "isolation with migration".

Three files were required as input: the joint site frequency spectra (sfs) for all six combinations of four populations, a template file () and a estimation file 
