# A counterexample to a conjecture of Plotnick on knot exteriors

This repository contains code and data accompanying the forthcoming preprint
*A counterexample to a conjecture of Plotnick on knot exteriors* by Tommy Hofmann
and John Nicholson.

The computations construct elements $u,v\in\mathbb ZG$, where
$G=\mathrm{SL}_2(\mathbb F_5)$, satisfying

$$
\varepsilon(u)=\varepsilon(v)=49,\qquad
tu=ut,\qquad tv=vt,\qquad uv=vu=1+20N,
$$

with $t=\left(\begin{smallmatrix}1&1\\0&1\end{smallmatrix}\right)$ and
$N=\sum_{g\in G}g$. Thus the images of $u$ and $v$ are mutually inverse units in
$\mathbb ZG/\mathbb ZN$. This is the algebraic construction used to resolve
Plotnick's conjecture concerning the homotopy types of two knot exteriors.

## Reproducing the computations

The computations were tested with Julia 1.13.1, Oscar 1.8.2, and Hecke 0.39.22.
The exact dependency versions are recorded in
[`notebooks/Manifest.toml`](notebooks/Manifest.toml).

The script
[`notebooks/discover_counterexample_oscar.jl`](notebooks/discover_counterexample_oscar.jl)
implements the search described in the manuscript's section *Finding the unit
algorithmically*. Run it from the repository root with

```sh
julia --project=notebooks notebooks/discover_counterexample_oscar.jl
```

An interactive version is available as
[`notebooks/discover_counterexample_oscar.ipynb`](notebooks/discover_counterexample_oscar.ipynb).
Open it with [IJulia](https://julialang.github.io/IJulia.jl/stable/) and run the
cells in order.
