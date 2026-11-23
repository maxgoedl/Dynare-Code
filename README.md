# Dynare-Code

A collection of Dynare mod-files for macro models. The goal is to provide transparent, self-contained replication files that document every modelling and calibration choice and can serve as a starting point for further model development.

## Replicability notes

Each mod-file has a header that states which results of the original paper are replicated. It also records where and why the implementation differs from the paper (e.g. continuous vs. discrete time, solution method) and any errors or inconsistencies found in the original.

## Compatibility

The mod-files have been tested with Dynare 7.0. Earlier Dynare versions may not support some features used here (e.g. `shocks(overwrite)` or macro-processor loops).

## Acknowledgements

Parts of the code and documentation in this repository were developed with the assistance of Claude (Anthropic). All model implementations were checked against the original papers by the author, who is responsible for any remaining errors.

## Content

### AG_unbalanced_growth.mod

Discrete-time version of Acemoglu, Daron and Guerrieri, Veronica (2008): "Capital Deepening and Nonbalanced Economic Growth", *Journal of Political Economy*, 116(3), pp. 467-498 (NBER Working Paper 12475, 2006/2007).

Replicates the benchmark calibration of Section 4: Table 1 (relative output, relative employment and the capital share, 1948 vs. 2004) and Figure 3 (λ, κ, r and σ_K over 3000 years for ε = 0.4, 0.5, 0.6). The 1948 values in Table 1 are matched exactly. The changes up to 2004 come out at about two thirds of the paper's. The interest rate path differs from Figure 3, most likely because of the polynomial approximation used in the paper (see the header for details).

This mod-file shows how to:
- simulate a model **without a balanced growth path** with the perfect-foresight solver, by normalising with the trend of the asymptotically dominant sector and putting the remaining non-stationary driver into an exogenous variable
- use a **conditional steady state** as an approximate terminal condition when the economy only converges asymptotically
- set **initial conditions from data** (1948 NIPA values) by backing out the technology levels

## License

Some files build on templates from Johannes Pfeifer's [DSGE_mod](https://github.com/JohannesPfeifer/DSGE_mod) repository, which is released under the GNU General Public License v3. The mod-files in this repository are distributed under the same license.
