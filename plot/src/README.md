# Plot

Scripts and data for generating all manuscript figures.

## Structure

- **[data/](data)** - Input data for figure generation (benchmarking results are in the `resources/` directory at repository root)
- **[src/](src)** - R and shell scripts for main, extended, and supplementary figures

## Usage

Each figure subdirectory in `src/` follows a consistent structure:
- `.sh` script - Executes the R plotting script
- `.r` script - Plotting code
- Output subdirectory - Generated figures (PDF, PNG, EPS, SVG) and intermediate data files

### Running scripts

Execute shell scripts from their respective directories:

```bash
cd src/fig1_main
sh panel_c.sh
