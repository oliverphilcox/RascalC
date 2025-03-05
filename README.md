# RascalC
**A Rapid Sampler For Large Cross-Covariance Matrices in C++**

C++ code to simulate correlation function covariance matrices from large surveys, using a grid and jackknife based approach. This can be used to find covariances of (a) the angularly binned anisotropic 2PCF, (b) the Legendre-binned anisotropic 2PCF and (c) the Legendre-binned isotropic 3PCF in arbitrary survey geometries. For (a) we can also compute the jackknife covariance matrix, which can be used to fit our non-Gaussianity model. There is additionally functionality to compute multi-tracer cross-covariances for the 2PCF.

For full usage, see the [ReadTheDocs documentation](https://rascalc.readthedocs.io/en/latest).

Any usage of this code should cite [Philcox et al 2020](https://arxiv.org/abs/1904.11070) (for the angularly binned 2PCF) and [Philcox & Eisenstein 2019](https://arxiv.org/abs/1910.04764) (for the original/accumulated Legendre-binned 2PCF, and 3PCF).
We also encourage the code users to cite [Rashkovetskyi et al 2023](https://arxiv.org/abs/2306.06320) (describing the post-reconstruction methodology, interfacing with [DESI](https://desi.lbl.gov), updated covariance comparison and validation techniques) and [Rashkovetskyi et al 2025](https://arxiv.org/abs/2404.03007) (for new/projected Legendre-binned 2PCF, extended validation with DESI DR1 mocks and an overview of the method).

***New for version 2***: Legendre moment covariances and the 3PCF

***New for version 3***: Python interface and library/package ([installation](https://rascalc.readthedocs.io/en/latest/library/installation.html), [basic usage](https://rascalc.readthedocs.io/en/latest/library/essential-functions.html), [practical examples](https://rascalc.readthedocs.io/en/latest/library/tutorials-examples.html))

## Authors

- Oliver Philcox (Columbia / Simons Foundation)
- Daniel Eisenstein (Harvard)
- Ross O'Connell (Pittsburgh)
- Alexander Wiegand (Garching)
- Misha Rashkovetskyi (Harvard)

We thank Yuting Wang and Ryuichiro Hada for pointing out and fixing a number of issues with the code and its documentation.
We are particularly grateful to Uendert Andrade for finding a wide variety of improvements and bugs, and to Jiaxi Yu for feedback on the Python interface and Jupyter notebook tutorial!
