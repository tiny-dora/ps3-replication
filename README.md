# Replication of Binstok, Gould and Kaplan (2024): "Policies or Prejudices?"

Problem Set 3, Empirical Methods, University of Zurich, HS 2026.

**Paper:** Binstok, N., Gould, E. D. and Kaplan, T. (2024), "Policies or Prejudices? An Analysis of Antisemitic and Anti-Israel Views on Social Media and Social Surveys", *Journal of the European Economic Association* 22(5), 2329–2373, doi:10.1093/jeea/jvae013.

**Original replication package of the authors:** https://academic.oup.com/jeea/article/22/5/2329/7613578#supplementary-data
We re-implemented the ALLBUS part of the authors' Stata do-file (`German_Social_Survey_2016_dofile.do`) in R. The authors' code is not included in this repository.

## What is replicated
- **Replication** (ALLBUS 2016): Table 6 (means and SDs), Table 7 (OLS), Table 9 (third panel, all columns), Table 10 (IV, all columns, both outcomes) and Figure 7 of the paper. All coefficients and sample sizes are identical to the paper; standard errors match to the published precision (`output/paper_vs_ours.csv`, `output/table6_vs_paper.csv`, `output/table9_vs_paper.csv`).
- **Extension:** weak-instrument-robust inference for Table 10 (tF, Anderson–Rubin) and a sensitivity analysis for a direct effect of the instrument (plausibly exogenous instrument, Conley, Hansen and Rossi 2012).

## How to run
1. Download the ALLBUS 2016 data from GESIS (free registration): **ZA5252, version 1.0.0**, https://search.gesis.org/research_data/ZA5252?doi=10.4232/1.12837. Save the Stata file as `data/raw/ZA5252_v1-0-0.dta`. The GESIS licence does not allow redistribution, so the data are not included in this repository.
2. Open `ps3-replication.Rproj` in RStudio (this sets the working directory to the repository root; all paths are relative).
3. Install the packages once: `install.packages(c("haven", "sandwich", "ivreg", "ggplot2"))`
4. Run `code/00_run_all.R`. It runs the three scripts in order (01 → 02 → 03) and writes all tables and figures to `output/`.

| Script | What it does | Report |
|---|---|---|
| `code/01_clean_allbus.R` | Raw ALLBUS file → cleaned data (sample restrictions and variables as in the authors' do-file) | – |
| `code/02_tables.R` | Replication of Tables 6, 7, 9, 10 and Figure 7; comparison with the published values | Table 1 (`table_main.tex`), Figure 1 (`figure7.pdf`); Appendix B: Tables 3–6 (`table6.tex`, `table7_*.tex`, `table9.tex`, `table10_*.tex`) |
| `code/03_extension.R` | Extension: first stage, reduced form, tF and Anderson–Rubin intervals; sensitivity to a direct effect of the instrument | Table 2 (`ext_weakiv_unjust.tex`), Figure 2 (`ext_plausexog.pdf`) |

Conventions: standard errors as in the authors' Stata code (HC1 for OLS, HC0 for IV); each regression uses its own complete cases, so N differs across columns as in Stata; the GESIS microdata in `data/raw/` and `data/processed/` are not committed.

## Folder structure
```
code/            R scripts (run code/00_run_all.R)
data/raw/        ALLBUS 2016 raw data (ZA5252_v1-0-0.dta, see "How to run")
data/processed/  cleaned data, created by 01_clean_allbus.R
output/          tables (.tex) and figures (.pdf) of the report; comparison with the paper
                 (paper_vs_ours.csv, table6_vs_paper.csv, table9_vs_paper.csv) and the numbers
                 behind Table 2 and Figure 2 (ext_weakiv.csv, ext_plausexog.csv)
```

## Software
Tested on macOS with R 4.2.3; haven 2.5.3, sandwich 3.1-3, ivreg 0.6-8, ggplot2 4.0.4. `code/00_run_all.R` writes the full list to `output/sessionInfo.txt`. No random numbers are used.

tF critical values from Lee, McCrary, Moreira and Porter (2022), as tabulated in the R package `ivDiag` (Lal and Xu, MIT licence).
