# Rootfinders for complex analytic functions on a square

This repo contains the code in Fortran 77 accompanying the manuscript:

H. Zhang, and V. Rokhlin, "Finding roots of complex analytic functions via generalized colleague matrices." Advances in Computational Mathematics 50.4 (2024): 71.

It includes subroutines for finding all zeros of complex analytic functions on a square.
The main programs are `comrootsq.f` and `comrootadapsq.f` for the nonadaptive and adaptive versions respectively. The code can be run in extended precision mode.

Run the bash script `comrootsq` and `comrootadapsq` to execute the program.

The special polynomial basis for rootfinding is precomputed in extended precision and stored in `orth_poly_data_leg_fin.f`, together with the corresponding three-term recurrence coefficients. Moreover, the QR decomposition of this basis is also precomputed and stored in `leastsq_data_square_leg_fin.f`.


