# Plot Source

R and shell scripts for generating manuscript figures.

## Structure

- **fig1_main/, fig2_main/** - Main manuscript figures
- **extended_figs/** - Extended Data figures
- **supplementary_figs/** - Supplementary Data figures

Each figure subdirectory contains:
- `.sh` script - Executes the R plotting script
- `.r` script - Plotting code and data processing
- Output subdirectory - Generated figures (e.g. PDF, PNG) and intermediate statistics used in the paper text

## Usage

Navigate to the figure directory and execute the shell script:

```bash
cd fig1_main
sh panel_c.sh
```

For main figures with multiple panels, run each panel script separately. 
Final composite figures (e.g., `figure1.pdf`, `figure2.pdf`) were assembled manually.

## Notes

- Main figures may include additional input files (e.g., `genes.panel_d.txt` in **fig2_main/**)
- Figure 8 (Supplementary Data) contains pre-generated panels; statistics are in `../data/time_mem_benchmark/`
- All scripts assume execution from within their respective directories

