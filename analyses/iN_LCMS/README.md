# LC-MS Experiments on induced Neurons

## Experiment Description
### LC-MS Experiments on iNs

#### Biphasic Extraction

iPSC-derived neurons were washed once with cold PBS (Fisher; Cat No. MT21040CM) and lifted off the plate with a cell scraper in 1 mL cold PBS. Cells were centrifuged at 2,000 × g for 5 min. PBS was removed, and cells were resuspended in 2 mL cold methanol for biphasic extraction. Chloroform (Sigma 1.02444; 4 mL cold) was added to each vial and mixed by vortexing for 1 min. Water (Sigma WX0001; 2 mL cold) was added and vortexed for 1 min. Vials were placed in 50 mL conical tubes and centrifuged at 3,000 rcf for 10 min for phase separation. The lower (chloroform) phase was collected (3 mL per sample) and transferred to new vials.

> **Note:** When samples were prepared by the Harvard Center for Mass Spectrometry, cell pellets arrived in 500 µL methanol, were vortexed, then transferred to 8 mL glass vials. Each received an additional 1.5 mL methanol and 4 mL chloroform, vortexed, and incubated in an ultrasound bath for 10 min. After adding 2 mL water and vortexing, phase separation was achieved by centrifugation at 800 rcf for 10 min at 4 °C. Upper aqueous phases went to metabolomics; lower chloroform phases went to lipidomics.

At least one blank (no cells) was processed alongside each extraction and run through the LC-MS pipeline.

#### Cell Pellet Preparation for LC-MS Lipidomics

Performed by the Harvard Center for Mass Spectrometry. Samples were dried under N₂ flow to \~1 mL, transferred to microcentrifuge tubes, and evaporated to dryness. Dried pellets were resuspended in chloroform (volume scaled to biomass; min. 60 µL), then split into two equal aliquots for positive and negative ionization analyses. For positive-only runs, samples were resuspended in 20–25 µL (min.) without splitting. After resuspension, samples were centrifuged (max speed, 10 min or 18,000 rcf for 20 min at 4 °C), and supernatants transferred into microinserts for LC-MS.

#### Cell Pellet Preparation for LC-MS Metabolomics

Also at the Harvard Center for Mass Spectrometry. Samples were dried under N₂ to \~1 mL, moved to microcentrifuge tubes, and evaporated to dryness. Dried pellets were resuspended in 50% acetonitrile/water (volume scaled to biomass; min. 20 µL). After centrifugation at max speed for 10 min, 12 µL or 15 µL of each supernatant (batch-dependent) was transferred into microinserts. Remaining supernatants were pooled to make batch-specific QC samples.

#### Media Preparation for LC-MS Metabolomics

Media (100 µL) was added to microcentrifuge tubes containing 1 mL methanol and held at –20 °C for 2 h. Samples were centrifuged at 18,000 rcf for 20 min at –9 °C; supernatants were transferred to new tubes and dried under N₂. Dried samples were resuspended in 50 µL 30% acetonitrile/water with 2 mM medronic acid, centrifuged again at 18,000 rcf for 20 min at 4 °C, and supernatants were moved into glass microinserts for LC-MS.

---

### LC-MS Lipidomics

Performed at Harvard C-MS (adapted from Miraldi et al. 2013). An Orbitrap Exactive Plus (Thermo) coupled to an Ultimate 3000 LC (Thermo) was used in top-5 data-dependent MS/MS mode, both positive and negative ionization. Separation was on a Biobond C4 column (4.6 × 50 mm, 5 µm; Dikma) at:

* **Flow**:

  * 0–5 min: 100 µL min⁻¹, 0% B
  * 5–55 min: ramp to 400 µL min⁻¹, linear gradient B 20→100%
  * 55–63 min: 500 µL min⁻¹, 100% B
  * 63–70 min: 500 µL min⁻¹, 0% B (re-equilibration)

* **Buffers (positive):**

  * A: 5 mM ammonium formate, 0.1% formic acid, 5% methanol in water
  * B: 5 mM ammonium formate, 0.1% formic acid, 5% water, 35% methanol in isopropanol

* **Buffers (negative):**

  * A: 0.03% ammonium hydroxide, 5% methanol in water
  * B: 0.03% ammonium hydroxide, 5% water, 35% methanol in isopropanol

Lipids were identified and integrated using LipidSearch v4.2.27, with manual curation before export.

---

### LC-MS Metabolomics

Also at Harvard C-MS, using a Vanquish LC + ID-X MS (Thermo). Samples (5 µL) were run on a ZIC-pHILIC PEEK-coated column (150 × 2.1 mm, 5 µm; Sigma) at 40 °C:

* **Mobile phases:**

  * A: 20 mM ammonium carbonate, 0.1% ammonium hydroxide in water
  * B: 97% acetonitrile in water

* **Gradient:**

  1. 0–30 s: flow ramp 0.05→0.15 mL min⁻¹ at 93% B
  2. 30 s–19 min: 93→40% B
  3. 19–28 min: 40→0% B
  4. 28–33 min: hold 0% B
  5. 33–36 min: 0→93% B
  6. 36–45 min: re-equilibrate at 93% B

Flow was 0.15 mL min⁻¹ (except initial ramp). Data were acquired in polarity-switching MS1 at 120,000 resolution (AGC 1×10⁵; m/z 65–1000). MS2/MS3 was run on pooled samples via AcquireX DeepScan (five reinjections each for positive and negative). A targeted metabolites standard mixture was run immediately after samples for quantitation.

## Analysis Description

#### LC-MS Lipidomics Data Analysis

Lipids were identified and their signals integrated using LipidSearch™ (v4.2.27, Mitsui Knowledge Industry, University of Tokyo). Integrations and peak quality were manually curated. Peak areas were background-corrected by subtracting three times the median peak area measured in blank samples; any resultant negative values were set to zero.

* **Between-cell-line comparisons**: Welch’s *t*-test (`scipy.stats.ttest_ind`, `equal_var=False`).
* **Within-cell-line (treatment) comparisons**: Student’s *t*-test (`scipy.stats.ttest_ind`, `equal_var=True`), assuming equal variance.

#### LC-MS Metabolomics Data Analysis

Data were processed in Compound Discoverer 3.2 (Thermo Fisher Scientific). Metabolite IDs were assigned by:

1. **Level 1**: MS²/MS³ spectral matching against a local mzVault library plus retention‐time confirmation with pure standards.
2. **Level 2**: Spectral matching via mzCloud.

All IDs were manually reviewed. Peak areas were background-corrected (subtracting three times the median blank), with negatives set to zero. Median-centered peak areas were then scaled using `StandardScaler()` from scikit-learn before principal component analysis (PCA). Three samples flagged by the Harvard C-MS as having low overall intensities were excluded from downstream analyses.

#### Targeted LC-MS Metabolite Analysis in Media Samples

Peak areas of CDP, CDP-choline, and choline in media were compared. Solvent blanks confirmed that CDP and CDP-choline were undetectable, while choline appeared at levels several orders of magnitude below those in media samples.

## Data Availability
- The raw will be available through [Dryad]().
- The processed data used in the following analyses will be available through [Dryad]().
  
## Code Overview
- Run `./lipidomics/SUB15127_lipidomics_Y622.ipynb` for p.Tyr622* lipidomics-related plots.
- Run `./lipidomics/SUB12877_G2_original.ipynb` for the p.Glu50fs*3 lipidomics-related plots.
- Run `./lipidomics/SUB14737_lipidomics_choline.ipynb` for p.Tyr622* +/- Choline lipidomics-related plots.
- Run `./metabolomics/SUB15127_media_metabolomics_targeted.ipynb` for the targeted metabolomics experiments.
- Run `./metabolomics/SUB15127_cells_metabolomics_untargeted.ipynb` for the untargeted metabolomics on p.Tyr622* +/- Choline.
