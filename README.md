isofrac
=======

Tools for computing equilibrium Mg isotope fractionation factors from CASTEP phonon calculations and visualizing fractionation between phases.

These scripts and notebooks were used in the preparation of the manuscript *“Controls on Mg isotopic fractionation between deep mantle phases and relict signatures of a terrestrial magma ocean”* (Andrew M. Walker, Remco C. Hin, Tim Elliott); see the preprint at https://doi.org/10.31223/X5SJ1F. Releases of this code are archived on Zenodo at https://doi.org/10.5281/zenodo.16744738.CASTEP data files can be found at https://doi.org/10.5281/zenodo.16760721.

This code is from a research workflow and is not intended as polished, general-purpose software. Expect to adapt pieces manually and bring appropriate expertise.

## Requirements

- Python 3 with `numpy`, `scipy`, and `matplotlib`.
- CASTEP with the PHONONS utility available on your `PATH`. CASTEP is requires a valid license; see https://castep.org for licensing details.

## Typical workflow

1. Run a standard CASTEP DFPT calculation for the phase of interest to generate a `<seedname>.check` file containing the force constants matrix. This should be available along with the corresponding `<seedname>.cell` and `<seedname>.param`
files.
2. Use `castep_isotope_sub.py` to generate isotope-substituted cells, launch PHONONS, fit the lnβ(T) parameters, and optionally clean up intermediates:
   ```
   python castep_isotope_sub.py MgO --plot --cleanup
   ```
   The script reports the fitted A/B/C coefficients for `1000·lnβ = A/T⁶ + B/T⁴ + C/T²`, writes temperature tables, and generates graphs (saved if `matplotlib.use('Agg')` is triggered).
3. Compare phases or visualize Δ²⁶Mg between reference and target phases using `alpha_plotting.py`, which evaluates `castep_isotope_sub.ln_beta_function` for user-provided coefficients.

If you already possess the light/heavy isotope `.phonon` files (e.g., from a previous run), you can bypass the automated CASTEP run and call `calc_beta.py` directly:
```
python calc_beta.py phase__isotope_l.phonon phase__isotope_h.phonon
```
This prints β, 1000·lnβ, and vectorized arrays for a default temperature grid.

Helper scripts such as `bulk_run_phonons.py` and `process_PVT_castep.py` illustrate how to automate multiple CASTEP runs (e.g. calculate or β for multiple pressures) and process equation-of-state data (needed to incorporate thermal expansion).

## Jupyter notebooks

These generate the figures in the manuscript.

- `alpha_beta_plots.ipynb` – interactive plotting of Δ²⁶Mg(T) using stored lnβ parameters and the helper plotting routines. Generates Figures 1, 2 and 3.
- `potential_figure.ipynb` – derives and plots the ionic-model generating Figures 4 and 8.
- `alpha_geotherm.ipynb` – full workflow from DFPT/structure inputs to Br–melt Mg fractionation along the liquidus. Generates Figures 5, 6, 7 and 9.
- `magma_ocean_model.ipynb` – models magma-ocean evolution scenarios and tracks associated isotopic signatures. Generates Figure 10

Each notebook assumes a Jupyter environment with `%matplotlib inline` enabled and can be executed independently once the required CASTEP outputs or fitted parameters are available.

## License

Distributed under the 3-clause BSD license; see `LICENSE` for the complete terms.
