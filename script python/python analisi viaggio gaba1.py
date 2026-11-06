import MDAnalysis as mda
from MDAnalysis.analysis import distances
import numpy as np
import matplotlib.pyplot as plt
import os

# 1. Caricamento topologia e traiettoria
topo = "step5_input.pdb"
traj = "step7_production.xtc"

print(f"Caricamento del sistema: {topo} e {traj}...")
u = mda.Universe(topo, traj)

# 2. Selezione dei due ligandi GABA
gaba_all = u.select_atoms("resname GAB or resname GABA")
resids = np.unique(gaba_all.resids)
print(f"Residui GABA identificati: {resids}")

if len(resids) < 2:
    raise ValueError("Attenzione: identificati meno di due residui GABA. Verifica il resname o resid nel PDB.")

# Assegnazione dei due GABA
gaba1 = u.select_atoms(f"(resname GAB or resname GABA) and resid {resids[0]}")
gaba2 = u.select_atoms(f"(resname GAB or resname GABA) and resid {resids[1]}")
protein = u.select_atoms("protein")

# Memorizziamo la posizione iniziale del Sito 1 e Sito 2 (frame 0) come array 2D (1, 3)
u.trajectory[0]
site1_ref_pos = np.array([gaba1.center_of_mass()])
site2_ref_pos = np.array([gaba2.center_of_mass()])

# 3. Analisi della traiettoria frame per frame
n_frames = len(u.trajectory)
times_ns = np.zeros(n_frames)
dist_from_site1 = np.zeros(n_frames)
dist_from_site2 = np.zeros(n_frames)
dist_gaba1_gaba2 = np.zeros(n_frames)
contacts_gaba1_protein = np.zeros(n_frames)

print(f"Analisi cinematica su {n_frames} frame con correzione PBC...")

for i, ts in enumerate(u.trajectory):
    times_ns[i] = ts.time / 1000.0  # tempo in ns
    
    # Coordinate COM come array 2D (1, 3) per distance_array con box PBC
    com_g1 = np.array([gaba1.center_of_mass()])
    com_g2 = np.array([gaba2.center_of_mass()])
    
    # Distanze corrette per PBC dal riferimento iniziale delle tasche
    dist_from_site1[i] = distances.distance_array(com_g1, site1_ref_pos, box=ts.dimensions)[0, 0]
    dist_from_site2[i] = distances.distance_array(com_g1, site2_ref_pos, box=ts.dimensions)[0, 0]
    
    # Distanza reciproca tra i due GABA corretta per PBC (minima immagine)
    dist_gaba1_gaba2[i] = distances.distance_array(com_g1, com_g2, box=ts.dimensions)[0, 0]
    
    # Contatti superficiali: atomi della proteina entro 4.5 Å da GABA-1
    d_mat = distances.distance_array(gaba1.positions, protein.positions, box=ts.dimensions)
    contacts_gaba1_protein[i] = np.sum(np.min(d_mat, axis=1) <= 4.5)

# 4. Generazione del grafico
fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(10, 8), sharex=True, dpi=300)

# Pannello 1: Allontanamento da Tasca 1 e avvicinamento a Tasca 2
ax1.plot(times_ns, dist_from_site1, color="#1f77b4", lw=1.8, label="Distanza GABA-1 da Sito 1")
ax1.plot(times_ns, dist_from_site2, color="#ff7f0e", lw=1.8, label="Distanza GABA-1 da Sito 2")
ax1.axvline(7.8, color="gray", ls="--", alpha=0.7, label="Distacco iniziale (~frame 283)")
ax1.axvline(52.6, color="purple", ls="--", alpha=0.7, label="Collisione/Bussata (~frame 1900)")
ax1.set_ylabel("Distanza da Tasca (Å)", fontsize=10, fontweight="bold")
ax1.set_title("Cinematica del GABA-1: Distacco e Trasferimento Inter-Tasca (Corretto PBC)", fontsize=12, fontweight="bold")
ax1.grid(True, ls=":", alpha=0.6)
ax1.legend(loc="upper left", fontsize=9)

# Pannello 2: Distanza reciproca GABA-1 / GABA-2
ax2.plot(times_ns, dist_gaba1_gaba2, color="#d62728", lw=2.0)
ax2.set_ylabel("Distanza GABA1-GABA2 (Å)", fontsize=10, fontweight="bold")
ax2.set_title("Distanza Relativa tra i Due Neurotrasmettitori (Clash Sterico)", fontsize=11, fontweight="bold")
ax2.grid(True, ls=":", alpha=0.6)

# Pannello 3: Contatti superficiali con la proteina
ax3.plot(times_ns, contacts_gaba1_protein, color="#2ca02c", lw=1.5)
ax3.set_ylabel("Atomi Proteici a Contatto", fontsize=10, fontweight="bold")
ax3.set_xlabel("Tempo (ns)", fontsize=11, fontweight="bold")
ax3.set_title("Adesione Superficiale all'ECD (Electrostatic Channeling)", fontsize=11, fontweight="bold")
ax3.grid(True, ls=":", alpha=0.6)

plt.tight_layout()

# Salvataggio con rimozione preventiva per evitare Errno 22 (file lock su Windows/Dropbox)
output_img = "crosstalk_gaba_kinetics_pbc.png"
if os.path.exists(output_img):
    try:
        os.remove(output_img)
    except Exception:
        pass

plt.savefig(output_img, dpi=300)
print(f"Analisi completata con successo! Grafico salvato in '{output_img}'.")