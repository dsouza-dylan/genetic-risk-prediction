# AI Use
The author primarily used Anthropic's Claude Code (Opus 4.7) as a pair-programming collaborator. The author would prompt specific tasks, after which Claude would propose code and explanations. Before committing any of these, the author would review, adapt, and test every iteration locally. All analytical choices and verification of results were made by the author.

## How Claude was used

Prompts clustered into 5 categories across the project.

### 1. Repository and notebook structure
Folder layout, naming conventions, consistent section headers across notebooks.

Representative prompts:
- *"Suggest a repository layout for cis-eQTL discovery + GWAS PRS"*
- *"Give me a header for the data wrangling notebook using this format: ..."*
- *"Based on the data wrangling notebook, how should I structure..."*

### 2. Code generation and refactoring
Helper functions (vectorized OLS, expression PRS cross-validation with null model, GWAS PRS pipeline). Consolidation of 7 modules into one `helpers.py`. Removal of unused helpers.

Representative prompts:
- *"Port my vectorized OLS in R to Python"*
- *"Write a function that implements 60/30/10 CV for expression PRS with a null model"*
- *"I just cleaned the helpers.py file. Are there any helpers that are no longer used in the notebooks?"*

### 3. Biological and pedagogical writeups
Markdown cells explaining Geuvadis, cis-eQTLs, clumping + thresholding, MAF, PRS computation at an introductory level.

Representative prompts:
- *"Walk me through writing up Geuvadis expression data for NBK 1 based on this linked tutorial..."*
- *"Simplify this by explaining clumping vs thresholding in plain English"*
- *"Dumb down the English in these 4 steps"*

### 4. Debugging, sanity-checking, and auditing
Verifying outputs against expected numbers and catching bugs.

Representative prompts:
- *"Is this expected: Running genome-wide cis-eQTL scan..."*
- *"Does this graph check out given the analysis in NBK 4"*
- *"Why did my genome-wide results run so quickly in Python?"*

### 5. Writing guidance
README drafting, AI contribution statement, prose refinement.

Representative prompts:
- *"Give me a 2-3 line summary on my code that I can use for the README"*
- *"How should I structure an AI use markdown file?"*

## Key decision points

Specific prompts that materially shaped the repository's final form:

| Prompt | Outcome |
|---|---|
| *"Can we minimize helpers.py by saving intermediates?"* | Final architecture: a single `helpers.py` with no intermediate files, notebooks as the single orchestration layer |
| *"Can I add the genome-wide loop to NBK 2 itself"* | Cached genome-wide scan added at the end of NBK 2 |
| *"22 CSV files for genome-wide or one genome-wide CSV file?"* | Per-chromosome output split (`output/eqtl_results/chr{N}_eqtl_results.csv`) |
| *"Should I commit the data for reproducibility, given the number of sources?"* | Data-first commit strategy for the repository |

## Testing approach

Every committed piece of code was verified through at least one of:

- **End-to-end notebook execution** after any structural or logic change.
- **Cross-consistency checks between outputs:** the Chromosome 22 in-memory scan in NBK 2 returned 162,845 SNP-gene tests, exactly matching the Chromosome 22 slice of the genome-wide scan (162,845 rows after splitting). Sum of the 22 per-chromosome output CSVs = 6,183,044 matched the reported genome-wide total.
- **Null-model validation for expression PRS:** real vs permuted effect sizes (e.g., FAM118A test R² = 0.58 vs null R² = 0.07 → real signal; several other genes beat null, some did not).

## Scope and limitations

- Claude's suggestions were not committed without critique. Failed suggestions were corrected on review. For instance, a "triple-read" bug where the same genotype file was read three times instead of unpacked once was caught during a code audit.
- All interpretation of results, gene selection, and methodological framing were made by the author.
- All architectural decisions (data commit strategy, per-chromosome splits, gene sets, PRS thresholds, individual to score) were made by the author.