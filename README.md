# RascalC
**A Rapid Sampler For Large Cross-Covariance Matrices in C++**

C++ code to simulate correlation function covariance matrices from large surveys, using a grid and jackknife based approach. This can be used to find covariances of (a) the angularly binned anisotropic 2PCF, (b) the Legendre-binned anisotropic 2PCF and (c) the Legendre-binned isotropic 3PCF in arbitrary survey geometries. For (a) we can also compute the jackknife covariance matrix, which can be used to fit our non-Gaussianity model. There is additionally functionality to compute multi-tracer cross-covarainces for the 2PCF. (Code developed by Alexander Wiegand, Daniel Eisenstein, Ross O'Connell and Oliver Philcox)

For full usage, see the ReadTheDocs [documentation](https://rascalc.readthedocs.io/en/latest).

Any usage of this code should cite [Philcox et al. 2019](https://arxiv.org/abs/1904.11070) (for the angularly binned 2PCF) and [Philcox & Eisenstein 2019](https://arxiv.org/abs/1910.04764) (for the Legendre-binned 2PCF and 3PCF).

***New for version 2***: Legendre moment covariances and the 3PCF
