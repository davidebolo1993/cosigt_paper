# Plot

Scripts and data for generating all manuscript figures.

## Structure

- **[data/](data)** - Input data for figure generation (e.g., benchmarking results)
- **[src/](src)** - R and shell scripts for main figures, extended figures, and supplementary figures

## Usage

Each figure subdirectory in `src/` contains:
- `.sh` script - Generates the figure by calling the R script
- `.r` script - Contains plotting code
- Output subdirectory with figures in multiple formats (e.g. PDF, PNG) and intermediate data files

### Running scripts

Execute shell scripts from their respective directories:

```bash
cd src/fig1_main
sh panel_c.sh
```