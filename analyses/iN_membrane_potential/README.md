# Imaging Experiments on induced Neurons

## Experiment Description

## Data availability
- Raw data are available on [OSF](https://osf.io/htb23/).
- Processed data are available on [OSF](https://osf.io/mnysb/files/osfstorage).

## Code overview
- Run `tmrm_011825_011925_with_fccp.ipynb` for TMRM analysis + FCCP
- Run `tmrm_012625_with_choline.ipynb` for snaps of all conditions
- Run `combined_quantifications_MitoHealth_Vehicle_updated_plots` to plot MitoHealth intensities at baseline [^1]
- Run `combined_quantifications_MitoHealth_Treatment_updated_plots` to plot MitoHealth intensities after treatment[^1]

[^1] For fixed z-stack images, NeuN-positive cell bodies were segmented in 3D using the pre-trained "cyto2" model (Cellpose). Segmentation quality was manually verified (blinded), and low-quality images were excluded. Cell-level fluorescence intensities were computed as probability-weighted sums of voxel intensities, using segmentation-derived voxel probabilities. The code for this quantification can be found at [GitHub](https://github.com/djunamay/confocalQuant/tree/main). 
