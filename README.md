# CEM_dissertation

The code in this package was developed as apart of my dissertation at The University of Manchester entitled "Shrinking For Restoring Definiteness" which is set to appear on MIMS EPrint. The abstract can be found at http://wp.me/p4w2Pn-5D

Some code was not written by myself and I have included the relevant licenses in the source code of those functions.

# Reproducing Results

The experiments in my dissertation were performed using MATLAB R2015a running Fedora 21 Linux. A single thread per physical core was used and thus no hyper-threading. Reproducing results may be difficult in future releases of MATLAB since we use the function maxNumCompThreads to set the maximum number of computational threads in the experiments which is to be removed.

# Acknowledgements

I'd like to thank Rémi Bazin for his help with GitHub.
