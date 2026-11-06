# Genetic Risk Prediction

Four notebooks on the 1000 Genomes EUR reference panel: 

**Part 1** maps cis-eQTLs (±500 kb) in Geuvadis LCL expression and builds an expression PRS evaluated with cross-validation against a null model.

**Part 2** builds a GWAS-based polygenic risk score for height using clumping + thresholding at seven p-value cutoffs.


## Setup

Requires **Python 3.11+** and **plink2**.

```bash
pip install -r requirements.txt
```

Download plink2 from <https://www.cog-genomics.org/plink/2.0/> and place it at `tools/plink2` (or `tools/plink2.exe` on Windows).


## Run

```bash
cd analysis
jupyter lab
```

Run the notebooks in order:
`data_wrangling.ipynb` → `cis_eqtl_discovery.ipynb` → `expression_prs.ipynb` → `gwas_prs.ipynb`


## AI and human contribution
Claude Code (Opus 4.7) was used to improve the code readability, notebook structure, as well as for drafting helper functions. These outputs were human-verified at every iteration, and the committed pieces of code reflect a byproduct of pair-programming with Claude Code, wherein the author set direction through prompts and iteratively refined Claude's suggestions through review and in-notebook testing. All analytical choices, such as gene selection, methodological decisions, and repository architecture were made by the author.