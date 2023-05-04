# nufft1d_optimization

We try to find a faster algorithm for Non-uniform Discret Fourier Transform (NUDFT) by solving the corresponding optimization problem. 
And we compare its performance with the existing two algorithms, which are based on Guassian interpolation and low rank approximation respectively.
Codes are focus on one-dimension type-II NUDFT.

## Experiments

Main results are displayed in **demo.mlx**.

## Contact

Any discussion and improvements are very welcome! (yongqijin0112 at gmail.com)

## References

[1] Leslie Greengard and June Yub Lee. Accelerating the nonuniform fast Fourier transform. SIAM Review, 46(3):443–4.
https://doi.org/10.1137/S003614450343200X

[2] Diego Ruiz-Antolin and Alex Townsend. A nonuniform fast Fourier transform based on low rank approximation. jan 2017.
https://doi.org/10.1137/17M1134822

[3] nufft1d codes: https://github.com/davidkrantz/nufft1d
