# Toroidal Density-Equalizing Map

<img src = "https://github.com/garyptchoi/toroidal-density-equalizing-map/blob/main/cover.jpg" height="300" />

* Toroidal density-equalizing map (TDEM): Compute a density-equalizing map on a prescribed toroidal surface.

Any comments and suggestions are welcome. 

If you use this code in your work, please cite the following paper:

S. Yao and G. P. T. Choi,
"[Toroidal Density-Equalizing Map for Genus-One Surfaces.](https://arxiv.org/abs/2410.16833)"
Preprint, arXiv:2410.16833, 2024.

Copyright (c) 2024, Shunyu Yao, Gary P. T. Choi

https://github.com/garyptchoi/toroidal-density-equalizing-map

===============================================================

TDEM usage:
* `map = TDEM(v,f,population,R,r,dt,epsilon,max_iter)`

Input:
* `v`: nv x 3 vertex coordinates of a sliced toroidal surface
* `f`: nf x 3 triangulations of a sliced toroidal surface
* `population`: nf x 1 positive quantity
* `R`: the major radius of the torus
* `r`: the minor radius of the torus
* `dt`: step size (default = 0.1)
* `epsilon`: stopping parameter (default = 1e-3)
* `max_iter`: maximum number of iterations (default = 1000)

Output:
* `map`: nv x 3 vertex coordinates of the toroidal density-equalizing map
* `uv`: nv x 2 vertex coordinates of the corresponding planar map
