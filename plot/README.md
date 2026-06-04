# Plot

Scripts, input tables, and generated outputs for the manuscript figures.

## Structure

- **[data/](data)** - Figure-specific input tables that are not already under `../resources/`
- **[src/](src)** - R and shell scripts for main and supplementary figures, plus committed generated outputs

## Usage

Each figure subdirectory in `src/` contains:

- `.sh` script - Calls the matching R script with the intended input tables
- `.r` script - Performs plotting and any panel-level summary calculations
- Output subdirectory - Generated figures in PDF, SVG, PNG, and/or EPS plus intermediate statistics

### Running scripts

Execute shell scripts from their respective figure directories:

```bash
cd src/fig1_main
sh panel_c.sh
```

The scripts assume their current working directory is the directory containing
the `.sh` file. Main composite figures such as `figure1.pdf` and
`figure2_new.pdf` were assembled manually from generated panel outputs.
