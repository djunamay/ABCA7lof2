
## SUB12877

### Lipidomics Experiment Summary
[Summarized with ChatGPT model 40 (GPT-4o) based on original descriptions [here](https://storage.googleapis.com/abca7lof/LCMS/SUB12877/5341.SUB12877_Lipidomics.pdf)]

The lipidomics experiment was conducted using a **ThermoFisher QE+ mass spectrometer** coupled with an **Ultimate 3000** system. The experiment followed two methods: **Lipidomics_top5_NEGATIVE** and **Lipidomics_top5_POSITIVE**, using **HESI (+ or -) source**. The MS parameters included **Full MS to ddMS2 scan from 0 to 60 min**, with a **70k resolution, AGC target of 1e6**, and an **m/z range of 150 to 2000**. The **dd-MS2 settings** were **35k resolution, AGC target of 1e5, loop count of 5, and stepped NCE of 20, 30, 40**. Waste was diverted during **0-1.5 min and 60-70 min**. The **liquid chromatography (LC) method** utilized a **BioBond 5 µm C4 column (50 × 4.6 mm)** with a **column temperature of 25°C** and **sample temperature of 4°C**. Samples were injected at **20 µL** using a gradient of **methanol, water, and isopropanol** with solvent modifications for positive and negative modes. 

Samples were **prepared by the client**, dried under **nitrogen (N₂) flow**, and resuspended in **chloroform**. The lowest biomass sample was resuspended in **60 µL**, and other samples were scaled accordingly. Samples were **split into two aliquots**, one for each polarity.

### Data Information

| Data Type          | Description | Link |
|--------------------|-------------|------|
| **Raw Data**      | Raw LC-MS data (SUB12877) | [Download](https://www.dropbox.com/scl/fo/zwxx32par2phdua0pwp2z/AIT-bmbxfCj1h-EMuCpUjY4?rlkey=gk5qa1hybb7ckz6gyqwv5fdeq&st=qkm7ba8q&dl=0) |
| **Lipid Data CSV**  | Processed lipidomics dataset | [Download](https://storage.googleapis.com/abca7lof/LCMS/SUB12877/1096.SUB12877_lipidXData.csv) |
| **Lipid Info Excel** | Lipidomics metadata and annotations | [Download](https://storage.googleapis.com/abca7lof/LCMS/SUB12877/3497.SUB12877_LipidomicsData.xlsx) |

### Sample Descriptions

| Sample ID | Description |
|-----------|------------|
| **Blank Samples** | 2 Blank samples |
| **S1** | Control, No Treatment |
| **S2** | ABCA7 LoF (p.Glu50fs*3), No Treatment |
| **S3** | Control, H₂O |
| **S4** | ABCA7 LoF (p.Glu50fs*3), H₂O |
| **S5** | ABCA7 LoF (p.Glu50fs*3), Choline |

### Code for plotting 
`SUB12877.ipynb`

## SUB13853

### Metabolomics Experiment Summary
[Summarized with ChatGPT model 40 (GPT-4o) based on original descriptions [here](https://storage.googleapis.com/abca7lof/LCMS/SUB13853/7177.SUB13853__metabolomics.pdf)]

Samples were prepared by drying under nitrogen (N₂) flow, transferring to microcentrifuge tubes, and evaporating to dryness. They were then resuspended in **50% acetonitrile in water** (volumes scaled to biomass), centrifuged at **max speed for 10 minutes**, and **15 µL of supernatant** was transferred to microinserts, with the remaining volume pooled for a QC sample. The experiment was conducted on a **Thermo ID-X Tribrid mass spectrometer (Vanquish LC system)** using the **SaltyHILIC_A2B2_150ul_100C_65to1k_45min_MS1_pn** method. MS parameters included **HESI (+/-) source, 120k resolution, RF lens 30%, AGC target 25%, max IT 50ms, and m/z range of 65-1000**, with **internal calibration** and waste diversion at **0-1 min and 32-45 min**. LC analysis was performed using a **Zic pHILIC column (150 × 2.1 mm, 5 µm)** at **40°C column temperature and 4°C sample temperature**, with an **injection volume of 5 µL**. The mobile phases included **20 mM ammonium carbonate with 0.1% ammonium hydroxide in water (A) and 97% acetonitrile in water (B)**, with a gradient from **93% B to 0% B and back to 93%** over **45 minutes**. Samples were **normalized by cell count and median-centered in Compound Discoverer (CD)**. The LC-MS strategy included **no internal standards (IS), a sequential run order, no QC in small batches, and an AquireX deep scan MS2/MS3 approach with 5 levels in positive and negative mode**.

Compound Discoverer (CD) **version 3.3** was used for data processing, applying a series of analytical steps. **Feature extraction** identified peaks from MS1 data, grouping adducts into compounds. **Retention time alignment** ensured consistency across files, while **gap filling** integrated missing peaks if detected at low levels. **Background subtraction** removed peaks found in blanks. **Normalization** was performed using **SERRF (for large datasets with QC samples)** or **median centering** to correct systematic errors. **Compound identification** involved MS/MS matching against **mzCloud, mzVault, and ChemSpider** databases, with manual curation to validate results. Confidence levels were assigned, ranging from **Level 1 (high-confidence matches) to weaker candidates identified through predictive modeling**. **Data visualization and significance analysis** included log2 fold change, p-values, and statistical comparisons, with flagged features for background interference, poor integration, or ambiguous identification.

The blanks showed **high signals for many compounds**, leading **Compound Discoverer (CD)** to tag most as background. However, compounds with **library matches** were manually reviewed, and those with **confirmed IDs** were retained in the final dataset but flagged as **potential background**. **Normalization** using **cell number-based resuspension** was effective in balancing data across samples, ensuring a **similar distribution**. **Median centering** further refined the dataset by smoothing out minor remaining variations, improving overall comparability.

### Data Information

| Data Type          | Description | Link |
|--------------------|-------------|------|
| **Raw Data**      | Raw LC-MS data (SUB13853) | [Download]() |
| **Processed Metabolite Data**  | Processed metabolite data with "Targeted", "Untargeted", and "SampleNames" tabs| [Download](https://storage.googleapis.com/abca7lof/LCMS/SUB13853/1338.SUB13853_MetabolomicsData.xlsx) |

### Code for plotting 
`SUB13853.ipynb`

## SUB14737

### Lipidomics and Metabolomics Experiment Summary
- See above for summaries on lipidomics and metabolomics experiments (SUB12877 and SUB13853).
- Also see here for the original description of this experiment: [here](https://storage.googleapis.com/abca7lof/LCMS/SUB14737/1757.SUB14737__metabolomics_lipidomics.pdf)]
- Note that this experiment was run only in positive mode.

### Data Information

| Data Type          | Description | Link |
|--------------------|-------------|------|
| **Raw Data**      | Raw LC-MS data (SUB14737) | [Download]() |
| **Processed Metabolite Data**  | Processed metabolite data with "Targeted", "Untargeted", and "SampleNames" tabs| [Download](https://storage.googleapis.com/abca7lof/LCMS/SUB14737/7372.SUB14737_MetabolomicsData.xlsx) |
| **Lipid Data CSV**  | Processed lipidomics dataset | [Download](https://storage.googleapis.com/abca7lof/LCMS/SUB14737/5041.SUB14737_LipidXData.csv) |
| **Lipid Info Excel** | Lipidomics metadata and annotations | [Download](https://storage.googleapis.com/abca7lof/LCMS/SUB14737/4459.SUB14737_LididomicsData.xlsx) |


### Sample Descriptions

| Sample ID | Description |
|-----------|------------|
| **S1** | ABCA7 LoF (p.Tyr622*), H₂O |
| **S2** | ABCA7 LoF (p.Tyr622*), CDP-Choline |
| **S3** | Extraction Blank |

### Code for plotting 
`SUB14737.ipynb`

## Code overview
- Run `./plotting_inputs.ipynb`. to get the lipidomic input object.
- Run `./lipidomics_by_subclass.ipynb` to plot lipidomics aggregate data for the iPSC-neuron.
- Run `./figures.ipynb` to plot figures.
- Run `./metabolomics-iN.ipynb` to analyze iPSC-neuron metabolomics data.
- Run `./choline.ipynb` to analyze iPSC-neuron metabolomics data.
