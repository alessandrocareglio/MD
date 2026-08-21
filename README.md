## Repository Overview

This repository contains the data, scripts, and visualization sessions for the dynamical network and correlation analysis of the protein structure **PDB ID: 8DD2**.

### Files Included

*   **`Dynamic Network Analysis of 8DD2 in Complex with GABA and ATP.pse`**
    *Description:* A saved PyMOL session visualizing the protein scaffold. The residues are color-coded based on their normalized cross-correlation values (mapped onto the B-factor/tempfactor column) relative to the ATP-binding hub. It highlights the potential allosteric communication pathways and coupled movements induced by ligand binding.

*   **`network_analysis.py`**
    *Description:* Python script utilizing the `MDAnalysis` toolkit to perform trajectory structural alignment (on C$\alpha$ atoms) and compute the Dynamical Cross-Correlation Matrix (DCCM). It calculates the average correlation between the global protein residues and the specific core residues within 7.0 Å of the ATP binding site.

*   **`heatmap_generation.py`**
    *Description:* Python script using `NumPy` and `Matplotlib`/`Seaborn` to process the calculated DCCM data and generate a high-resolution correlation heatmap, visually representing the coupled motions and global dynamic network of the system.
