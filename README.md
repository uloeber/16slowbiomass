# Host Contamination in 16S rRNA Gene Sequencing

**Repository for the manuscript:**  
**"Effects of Host Contamination in 16S rRNA Gene Sequencing: Bias, Impact, and Mitigation Strategies"**

📌 **NCBI BioProject accession**: [PRJNA1268757](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1268757)  
📌 **Publication DOI (upon acceptance)**: *TBD*

---

## 📁 Repository Structure

```
├── input/                      # Metadata and processed feature tables
│   ├── METADATA.csv
│   ├── lotus2_SLV_raw/
│   ├── lotus2_SLV_cont_rm_raw/
│   ├── lotus2_SLV_raw_rmhumanlotus/
│   ├── lotus2SLV138_humancontUganda/
│   └── microbial_composition.txt
├── figures/                    # Output plots (svg/pdf)
├── ddqpcr/                     # Raw and processed ddPCR results
├── intermediate/              # Temporary or intermediate processed files
├── analysis.R                 # Main analysis script (all figures/stats)
├── ddpcr.R                    # Digital droplet PCR quantification visualization
├── humanContam.Rproj          # RStudio project file
└── README.md                  # This file
```

---

## 📊 Project Summary

Low-biomass samples are highly susceptible to off-target amplification of host DNA during 16S rRNA gene sequencing. This project:
- Simulates bacterial-to-human DNA ratios in mock communities
- Evaluates contamination via raw, pre-clustering, post-clustering, and `decontam` filtering
- Validates findings with digital droplet PCR (ddPCR)

🧪 Results show:
- Host DNA >90% severely distorts microbial profiles
- Total DNA concentration has limited biasing effect
- Post-clustering OTU filtering is most effective and computationally efficient
- `decontam` struggles with low-biomass data

---

## 🔧 Setup

### Prerequisites
Install the following R packages:

```r
install.packages(c("ggplot2", "vegan", "cowplot", "phyloseq", "decontam", "dplyr", "reshape2", "tidyr", "ggrepel", "ggpubr", "latex2exp", "stringr", "lme4", "lmtest"))
devtools::install_github("mikemc/metacal")
```

---

## 📈 Reproducing the Analysis

1. **Clone the repository:**

```bash
git clone https://github.com/uloeber/16slowbiomass.git
cd 16slowbiomass
```

2. **Open the R project:**
   Open `humanContam.Rproj` in RStudio.

3. **Run the main analysis:**
   Execute `analysis.R` to:
   - Read and clean data
   - Perform filtering comparisons (naive, pre/post-clustering, decontam)
   - Generate Figures 1–6 and supplementary figures

4. **Run ddPCR validation:**
   ```r
   source("ddpcr.R")
   ```

---

## 📂 Key Outputs

- `figures/`: Plots used in the manuscript
- `intermediate/percNoHit.tsv`: Table showing proportion of non-bacterial reads
- `ddqpcr/19124_all_ddqpcr_results.csv`: ddPCR quantification of 16S/18S rRNA genes

---

## 🧬 Data Availability

All raw sequencing data are publicly available via [NCBI SRA – PRJNA1268757](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1268757).  

Code and scripts in this repository are sufficient to fully reproduce the published results.

---

## 📝 Citation

Please cite the manuscript if you use this repository or its methods:

> Birkner T, Bartolomaeus TUP, McParland V, Forslund SK, Löber U. *Effects of Host Contamination in 16S rRNA Gene Sequencing: Bias, Impact, and Mitigation Strategies*. Microbiome (2025). [DOI pending]

---

## 👩‍🔬 Authors

- **Ulrike Löber** – [ulrike.loeber@mdc-berlin.de](mailto:ulrike.loeber@mdc-berlin.de)
- Sofia K. Forslund, Theda Bartolomaeus, Till Birkner, and others (see manuscript)

---

## 📜 License

MIT License. See `LICENSE` for details.

---
