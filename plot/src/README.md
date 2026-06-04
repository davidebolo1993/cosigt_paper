# Plot Source

R and shell scripts for generating manuscript figure panels.

## Structure

- **fig1_main/** - CYP2D6 graph and cluster-visualization panels
- **fig2_main/** - COSIGT benchmark summary panels
- **supplementary_figs/** - Expanded benchmark, HLA, and runtime/resource panels

Each figure subdirectory contains:

- `.sh` script - Executes the R plotting script
- `.r` script - Plotting code and data processing
- Output subdirectory - Generated figures and intermediate statistics used in the paper text

## Usage

Navigate to the figure directory and execute the shell script:

```bash
cd fig1_main
sh panel_c.sh
```

For main figures with multiple panels, run each panel script separately. Final
composite figures, such as `figure1.pdf` and `figure2_new.pdf`, were assembled
manually from panel outputs.

## Figure Map

- `fig1_main/panel_c.sh` and `panel_d.sh` plot the CYP2D6 region data in `../data/cyp2d6/`.
- `fig2_main/panel_b.sh` summarizes leave-zero-out COSIGT vs Locityper QV categories for CMRGs and SVs.
- `fig2_main/panel_c.sh` summarizes leave-all-out QV fraction relative to the best alternate-pangenome hit.
- `fig2_main/panel_d.sh` compares selected genes across modern and contaminated aDNA settings.
- `supplementary_figs/fig4_supp` through `fig8_supp` expand COSIGT vs Locityper per-gene comparisons across CMRG coverages and SVs.
- `supplementary_figs/fig9_supp` through `fig11_supp` expand leave-all-out QV-fraction comparisons.
- `supplementary_figs/fig12_supp` and `fig13_supp` summarize selected aDNA conditions and demographic simulation settings.
- `supplementary_figs/fig14_supp` plots mean COSIGT QV against region length.
- `supplementary_figs/fig15_supp` plots mean COSIGT QV against average normalized cluster distance.
- `supplementary_figs/fig18_supp` and `fig19_supp` plot HLA summary accuracy from `../data/hla/`.

## Notes

- Main figures may include additional input files, such as `genes.panel_d.txt` in `fig2_main/`.
- Some supplementary PDFs are committed as pre-generated panels without matching scripts in this directory.
- All scripts assume execution from within their respective directories.
