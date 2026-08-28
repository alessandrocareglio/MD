Markdown
# Positive Allosteric Modulation of the GABA_A Receptor by Extracellular ATP

This repository provides the simulation parameter files, analysis scripts, extracted trajectory frames, and 3D PyMOL visualization sessions accompanying the study:
> **"Positive Allosteric Modulation of the GABAA Receptor by Extracellular ATP: Structural Decoding, Network Dynamics, and Conformational Coupling"**

---
## 📂 Repository Contents

### **Visualization Sessions & Structural Models**
- **`8DD2_ATP_GABA_Allosteric_Trajectory_Session.pse`**  
  PyMOL multi-frame session displaying the MD trajectory of ATP binding and accommodation at the basal $\alpha_1/\gamma_2$ interface, along with the structural response of the proximal orthosteric GABA pocket.
- **`Dynamic Network Analysis of 8DD2 in Complex with GABA and ATP.pse`**  
  PyMOL session visualizing the 3D cryo-EM structure (PDB ID: 8DD2) color-coded by normalized cross-correlation values ($C_{ij}$) mapped onto the B-factor column, highlighting the continuous allosteric coupling pathway.
- **`8DD2_pathway_correlated.pdb`**  
  Coordinate file with computed DCCM allosteric correlation values stored in the temperature factor (B-factor) field for structural rendering.
- **`frames/`**  
  Directory containing discrete extracted trajectory PDB frames along the 1.0 ns production run.

### **Analysis & Workflow Scripts**
- **`dccm_analysis.py`**  
  Python script using `MDAnalysis` to perform backbone C$\alpha$ alignment and compute the Dynamic Cross-Correlation Matrix (DCCM) across all pentameric subunits.
- **`map_dccm_to_pdb.py`**  
  Script calculating the average correlation between all receptor residues and the 7.0 Å ATP-binding hub, writing the normalized scores into the PDB B-factor column.
- **`crea_video.py`**  
  Python utility using `imageio` to compile high-resolution rendering frames into smoothed, looping trajectory MP4 videos.

### **GROMACS / CHARMM-GUI Parameters**
- **`step5_production.mdp`**  
  GROMACS molecular dynamics production run parameter file (CHARMM36m force field, PME electrostatics, 303.15 K).
- **`topol.top`**  
  Complete system topology file generated via CHARMM-GUI.

---

## 🔬 Key Biophysical Findings
1. **Concerted Multi-Loop Capping:** The high negative charge density of ATP engages $\alpha_1$ Loops A–C and $\gamma_2$ Loops D'–F, acting as a concerted electrostatic sensor.
2. **Asymmetric Orthosteric Priming:** Tangential mechanical tension propagates across rigid $\beta$-sheets to relieve initial steric clashes at the proximal GABA site ($\beta_2^+/\alpha_1^-$), driving an induced-fit affinity transition ($\Delta G$ shifts from $+6.09$ to $-5.21\text{ kcal/mol}$).
3. **Pore Stabilization & Kinetic Damping:** Endoluminal projection of the triphosphate chain reorganizes pore electrostatics (net charge relaxes from $+6$ to $+3$) and dynamically damps M2 helices, providing a structural rationale for prolonged open-state dwell times ($\tau_{\text{open}}$).
