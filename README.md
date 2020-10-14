# MP2RAGE scripts

- T1 maps estimation using the MP2RAGE sequence
- T1 map correction using an additional B1 map
- Background noise removal by using a "robust"/regularized version of the combination of the two inversion time images

These functions can be called with BIDS-wrappers (the functions named `bids_*`) to automatically process BIDS data repositories
