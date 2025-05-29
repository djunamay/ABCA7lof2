# Seahorse Experiments on induced Neurons

## Experiment Description
iPSCs were differentiated as described above directly on 96-well Agilent Seahorse XFe96/XF Pro cell culture microplates and matured for 28 days before assaying on a Seahorse XFe96 Analyzer. Seahorse XF Cell Mito Stress Test and Oxidation Stress Tests were performed according to manufacturer protocol with the following final drug concentration: Oligomycin, 2.5 µM; FCCP, 1 µM; Rotenone/Antimycin, 0.5 µM. 

The oxygen consumption rate (OCR) of cells was determined over time using a Seahorse XF Analyzer. Prior to analysis, OCR curves were visually inspected in a blinded manner to exclude wells that did not respond to drug injections. To calculate per-well total oxygen consumption for a given experimental period (e.g., under basal conditions prior to injections of uncouplers), integrals between specific experimental time points were computed from the OCR curve. The following measurements were made:

1. **Basal respiration** was computed as the total oxygen consumption prior to oligomycin injection.
2. **Proton leak** was computed as the total oxygen consumed after oligomycin injection and prior to FCCP injection.
3. **Maximal respiration** was computed as the total oxygen consumption after FCCP and prior to Rotenone + Antimycin injection.
4. **Relative uncoupling** was computed as the fraction of basal respiration attributed to proton leak.
5. **Spare respiratory capacity** was determined as the ratio of maximal respiration to basal respiration.

## Data availability
- Raw and Processed data will be made available on Dryad

## Code overview
- Run `./seahorse.ipynb`.
