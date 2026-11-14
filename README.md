# Thermodynamic Asymmetry and Electromechanical Gating in the Human $\text{GABA}_\text{A}$ Receptor Primed by Extracellular ATP

This repository provides the computational workflows, Python analysis scripts, Jupyter notebooks, parameter files, and direct links to the full-length all-atom molecular dynamics (MD) trajectories accompanying the study:

> **"Thermodynamic Asymmetry and Electromechanical Gating in the Human $\text{GABA}_\text{A}$ Receptor Primed by Extracellular ATP"**  
> *Alessandro Careglio* — B.S. in Herbal Sciences and Technologies, PharmD Student, University of Turin (Università degli Studi di Torino)  
> ORCID: [0000-0002-3005-4689](https://orcid.org/0000-0002-3005-4689)

---

## 🌐 External Trajectory Files (100 ns Polarized POPC Bilayer Production Run)

Due to GitHub file-size limits, the full-length $100\text{ ns}$ GROMACS production trajectory ($5,001$ frames, explicit POPC lipid bilayer, constant electric field $V_m = -70\text{ mV}$) and starting system coordinates are hosted via persistent cloud links:

* **Starting System Coordinates & Topology Reference (`step5_input.pdb`):**  
  [Download `step5_input.pdb` via Dropbox](https://www.dropbox.com/scl/fi/dbxdaemts4q2kgj2jxesm/step5_input.pdb?rlkey=qtlxo1rf39sxe1tq5p7yay11r&st=vl6xksul&dl=0)
* **100 ns Production Trajectory (`step7_production.xtc`):**  
  [Download `step7_production.xtc` via Dropbox](https://www.dropbox.com/scl/fi/gyucrve8fmxvaiq1qv571/step7_production.xtc?rlkey=e26giijy6p9dmr28ou9x7jkzt&st=tygtqk2q&dl=0)

### Quick Visualization in PyMOL
Place both downloaded files in the same working directory and run the following commands in the PyMOL command line:
```python
load step5_input.pdb, gabaa_md
load_traj step7_production.xtc, gabaa_md, interval=5
hide everything
show cartoon, polymer.protein
show sticks, resn ATP+GABA

