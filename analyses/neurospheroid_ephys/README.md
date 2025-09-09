# Ephys on cortical organoids

## Experiment Description
Electrophysiological recordings were performed using an Axon Multiclamp 700B amplifier and Clampex 11.2 software (Molecular Devices). Cells were visualized using infrared differential interference contrast (IR-DIC) imaging (Olympus BX-50WI microscope), placed in a recording chamber, and perfused continuously at 2 mL/min ( 32°C) with oxygenated artificial cerebrospinal fluid (ACSF, containing 125 mM NaCl, 2.5 mM KCl, 1.2 mM NaH2PO4•H2O, 2.4 mM CaCl2•2H2O, 1.2 mM MgCl2•6H2O, 26 mM NaHCO3, and 11 mM D-Glucose).
Action potentials were elicited by injecting current steps in current-clamp mode. Whole-cell currents were recorded from a holding potential of -80 mV by stepping to various voltages in voltage-clamp mode. Spontaneous firing was recorded in cell-attached configuration. Recordings were filtered at 1 kHz (four-pole Bessel filter), digitized at 10 kHz with a Digidata 1550B interface (Molecular Devices). Pipette solution contained 120 mM K-gluconate, 5 mM KCl, 2 mM MgCl2•6H2O, 10 mM HEPES, 4 mM ATP, and 0.2 mM GTP. Data were analyzed using pClamp 11.2 and GraphPad Prism 10.
For electrophysiological recordings from cortical organoids, day 150 organoids were dissociated using Accutase (Stem Cell Technologies, #07920, 40 min, 37°C), plated onto #1 glass coverslips (Fisher Scientific, #50-194-4702) coated with PDL, laminin, and Matrigel, and maintained in 2D culture with or without 100 µM CDP-choline for two weeks prior to recordings.
Spontaneous action potential outliers were identified using the IQR method (values outside Q1-2×IQR or Q3+2×IQR) and removed, resulting in the exclusion of two data points (values: 9.38 in p.Tyr622*; 6.15 in p.Tyr622*+Choline). Cells recording zero spontaneous potentials (likely glial) were also excluded.

## Data availability
- Raw data are available on [OSF](https://osf.io/em92c/).

## Code overview
- Run `plotting.ipynb` to generate panels.
