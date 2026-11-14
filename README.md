# Thermodynamic Asymmetry and Electromechanical Gating in the Human GABA_A Receptor Primed by Extracellular ATP

This repository provides the computational workflows, Python analysis scripts, Jupyter notebooks, parameter files, and direct links to the full-length all-atom molecular dynamics (MD) trajectories accompanying the study:

> **"Thermodynamic Asymmetry and Electromechanical Gating in the Human GABAA Receptor Primed by Extracellular ATP"**  
> *Alessandro Careglio* — B.S. in Herbal Sciences and Technologies, PharmD Student, University of Turin (Università degli Studi di Torino)  
> ORCID: [0000-0002-3005-4689](https://orcid.org/0000-0002-3005-4689)

---

## External Trajectory Files (100 ns Polarized POPC Bilayer Production Run)

Due to GitHub file-size limits, the full-length 100 ns GROMACS production trajectory (5,001 frames, explicit POPC lipid bilayer, constant electric field Vm = -70 mV) and starting system coordinates are hosted via persistent cloud links:

* **Starting System Coordinates & Topology Reference (`step5_input.pdb`):**  
  [Download `step5_input.pdb` via Dropbox](https://www.dropbox.com/scl/fi/dbxdaemts4q2kgj2jxesm/step5_input.pdb?rlkey=qtlxo1rf39sxe1tq5p7yay11r&st=vl6xksul&dl=0)
* **100 ns Production Trajectory (`step7_production.xtc`):**  
  [Download `step7_production.xtc` via Dropbox](https://www.dropbox.com/scl/fi/gyucrve8fmxvaiq1qv571/step7_production.xtc?rlkey=e26giijy6p9dmr28ou9x7jkzt&st=tygtqk2q&dl=0)

### Quick Visualization in PyMOL
Place both downloaded files in the same working directory and run the following commands in the PyMOL command line:

    load step5_input.pdb, gabaa_md
    load_traj step7_production.xtc, gabaa_md, interval=5
    hide everything
    show cartoon, polymer.protein
    show sticks, resn ATP+GABA

*(Note: `interval=5` loads every 5th frame to optimize system memory while preserving smooth conformational dynamics; remove `interval=5` to load all 5,001 frames).*

---

## Repository Contents

### Primary Analysis Notebooks & Scripts
* **`MD_GABA_3.ipynb`**  
  Main reproducible Jupyter Notebook integrating trajectory post-processing, axial pore radius profilometry, spatiotemporal radial surface mapping R(z,t), luminal water density quantification, 3D vectorial Cl- tracking, side-chain rotameric dihedral profiling (chi1, chi2), DCCM computation, and `gmx_MMPBSA` energy parsing[cite: 64].
* **`script python/`**  
  Modular Python 3 scripts (built on `MDAnalysis`, `Bio3D`, `NumPy`, `Pandas`, and `Matplotlib`) for parsing `gmx_MMPBSA` decomposition matrices and generating publication-grade (300 DPI) figures[cite: 64].
* **`dccm_analysis.py`**  
  Python utility using `MDAnalysis` to perform backbone C-alpha alignment and compute the Dynamic Cross-Correlation Matrix (DCCM) across all pentameric subunits[cite: 64].
* **`map_dccm_to_pdb.py`**  
  Script mapping normalized cross-correlation scores (Cij) onto the PDB B-factor column for 3D structural rendering[cite: 64].

### Structural Models, PyMOL Sessions & Supplementary Movies
* **`Dynamic Network Analysis of 8DD2 in Complex with GABA and ATP.pse`**  
  PyMOL session visualizing the 3D cryo-EM structure (PDB ID: 8DD2) color-coded by normalized cross-correlation values (Cij) mapped onto the B-factor column, highlighting the continuous allosteric coupling pathway[cite: 64].
* **`8DD2_pathway_correlated.pdb`**  
  Coordinate file with computed DCCM allosteric correlation values stored in the temperature factor (B-factor) field[cite: 64].
* **`frames/`**  
  Directory containing discrete extracted trajectory PDB snapshots from preliminary isotropic relaxation runs[cite: 64].
* **Supplementary Trajectory Animations:**  
  * **Supplementary Movie 1** (Iris-like gating and pore bottleneck fluctuations): [Watch on YouTube](https://www.youtube.com/watch?v=HsiEkBDsSWs)[cite: 54]
  * **Supplementary Movie 2** (Conformational morphing along the gating coordinate): [Watch on YouTube](https://www.youtube.com/watch?v=MRronyy_uuc)[cite: 54]

### GROMACS / CHARMM-GUI Parameters
* **`step5_production.mdp`** — GROMACS production run parameter file (CHARMM36m force field, PME electrostatics)[cite: 64].
* **`topol.top`** — Complete system topology file generated via CHARMM-GUI[cite: 64].

---

## Key Biophysical Findings

1. **Kinetic Damping & Allosteric Priming:** Confinement of polyanionic ATP (-4e) at the extracellular beta2(+)/alpha1(-) interface—anchored primarily by **Lys156** (DeltaG_res = -14.2 kcal/mol)[cite: 31] and locally compensated by discrete **Asp156** rotameric switching (g- to trans, Delta_theta = 139.4 deg)[cite: 31]—suppresses interfacial loop fluctuations and triggers Loop C capping (13–16 ns)[cite: 31, 41].
2. **Thermodynamic Asymmetry Across Orthosteric Sites (MM/GBSA):** The primary orthosteric pocket (**GABA-1**, site 359) exhibits spontaneous, exergonic binding (DeltaG_bind = -2.54 kcal/mol; vacuum Coulombic attraction and polar desolvation cancel within +0.05 kcal/mol), led by **Arg67** (DeltaG_res = -3.08 kcal/mol) and the interfacial aromatic cage[cite: 31]. Conversely, the secondary pocket (**GABA-2**, site 360) populates a compressed, metastable state (DeltaG_bind = +39.08 kcal/mol) dominated by steric clashes (DeltaE_vdw = +33.42 kcal/mol)[cite: 31].
3. **ATP as a High-Elastic-Strain Biomechanical Wedge:** Allosteric ATP registers an apparent DeltaG_bind = +903.83 kcal/mol (DeltaE_GB = +639.76 kcal/mol; DeltaE_vdw = +496.59 kcal/mol), quantifying stored elastic strain transduced via the ECD beta-sandwich lever to drive a concerted pore dilation pulse at t = 18.80 ns[cite: 31, 52].
4. **Apical Cl- Recruitment & M2 Hydrophobic Vapor Lock:** 3D vectorial tracking of all 146 Cl- ions confirms axial recruitment through the central apical ECD funnel (15 entry events vs. 0 lateral fenestration events), while a continuous 38-A completely dewetted M2 lumen (r = 2.2 A, density = 0.0000 molecules/nm3) enforces the hydrophobic *vapor lock* characteristic of the closed/desensitized Cryo-EM state[cite: 31, 50].
