# Marine PETase Homolog Analysis (Tara Oceans)

Investigation of PETase homologs across global ocean metagenomic datasets from the Tara Oceans expedition.

## Repository Structure

```text
.
├── docs/
│   ├── commands.txt               # Raw MMseqs2 preprocessing commands
│   ├── main.pdf                   # Final compiled 8-page report
│   ├── main.typ                   # Typst report source
│   └── ref.bib                    # References
├── fig/                           # Generated publication figures
├── tools/                         # Analysis, statistics, and plotting scripts
├── etc/                           # Supplementary materials (ignored by git)
└── data/                          # Metagenomic datasets (ignored by git)
```

## Usage


### 1. Statistical Analysis & Calculations (Reproducibility)
Run all scripts from the repository root:

- **Reproduce Environmental Statistics & Correlations (Report Page 6):**
  ```bash
  python tools/calc_environmental_stats.py
  ```
  *Calculates depth layer distributions, dominant cluster parameters, Spearman rank correlations (temperature $\rho = -0.344$, depth $\rho = -0.183$, salinity $\rho = -0.250$), and vertical Kruskal-Wallis test ($H = 22.35, p = 1.40 \times 10^{-5}$).*

- **Reproduce 3D Structural Metrics & Permutations (Report Page 6–7):**
  ```bash
  python tools/calc_structural_stats.py
  ```
  *Calculates centroid distance correlation ($\rho = 0.33, p < 0.001$), $C_\beta$ 8-neighbor packing density correlation ($\rho = 0.50, p < 0.001$, 1,000 permutations), and solvent accessibility (RSA) contrasts on PDB 6QGC.*

- **Reproduce Sequence Filtering & Active Site Motifs (Report Page 3–4):**
  ```bash
  python tools/calc_sequence_stats.py
  ```
  *Calculates MMseqs2 hit filtering, identity range (19.5%–46.2%), reliable coverage positions (192 aa), G-x-S-x-G motif conservation, and catalytic triad invariance ($H = 0.00$ bits).*

### 2. Generating Figures
- **Filter search hits:**
  ```bash
  python tools/filter_search_out.py
  ```
- **Generate MSA and Shannon entropy plots** (Figures 2 & 5):
  ```bash
  python tools/reconstruct_msa_and_entropy.py
  ```
- **Generate dominant cluster maps & depth stratification** (Figures 3 & 4):
  ```bash
  python tools/plot_petase_clusters.py
  ```
- **Plot global oceanographic world map**:
  ```bash
  python tools/plot_world_map.py temperature   # options: temperature, salinity, depth, oxygen
  ```

### 3. Compiling the Report
Recompile the report from Typst source to PDF:
```bash
typst compile docs/main.typ docs/main.pdf
```
