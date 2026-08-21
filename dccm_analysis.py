import MDAnalysis as mda
from MDAnalysis.analysis import align
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# 1. Carica struttura e traiettoria
# Usa prod.tpr o prod.gro come topologia e prod.xtc per le coordinate
u = mda.Universe("prod.tpr", "prod.xtc")

# 2. Seleziona i soli atomi C-alfa della proteina
calphas = u.select_atoms("protein and name CA")
n_atoms = len(calphas)
print(f"Atomi C-alfa selezionati: {n_atoms}")

# 3. Allinea la traiettoria sul primo frame per rimuovere traslazioni/rotazioni
ref = u.select_atoms("protein and name CA")
aligner = align.AlignTraj(u, ref, select="protein and name CA", in_memory=True).run()

# 4. Estrai le posizioni dei C-alfa lungo tutti i frame
coords = []
for ts in u.trajectory:
    coords.append(calphas.positions)
coords = np.array(coords) # Dimensioni: (n_frames, n_atoms, 3)

# 5. Calcola le fluttuazioni medie (delta_r = r - <r>)
mean_pos = np.mean(coords, axis=0)
delta_r = coords - mean_pos # Dimensioni: (n_frames, n_atoms, 3)

# 6. Calcola la Matrice di Cross-Correlazione Dinamica (DCCM)
# C_ij = <delta_r_i . delta_r_j> / sqrt(<delta_r_i^2> * <delta_r_j^2>)
print("Calcolo della matrice DCCM in corso...")
dccm = np.zeros((n_atoms, n_atoms))

dot_prod_matrix = np.zeros((n_atoms, n_atoms))
for frame in range(delta_r.shape[0]):
    dr = delta_r[frame]
    dot_prod_matrix += np.dot(dr, dr.T)

mean_dot = dot_prod_matrix / delta_r.shape[0]
diag = np.diag(mean_dot)
dccm = mean_dot / np.sqrt(np.outer(diag, diag))

# 7. Genera e salva la mappa di calore (Heatmap)
plt.figure(figsize=(10, 8), dpi=300)
cmap = sns.diverging_palette(240, 10, as_cmap=True) # Da blu (-1) a rosso (+1)

sns.heatmap(
    dccm,
    cmap="coolwarm",
    vmin=-1.0,
    vmax=1.0,
    center=0.0,
    square=True,
    cbar_kws={'label': 'Cross-Correlation Coefficient (Cij)'}
)

plt.title("Dynamic Cross-Correlation Matrix (DCCM) - 8DD2", fontsize=14, pad=15)
plt.xlabel("C-alpha Residue Index", fontsize=12)
plt.ylabel("C-alpha Residue Index", fontsize=12)

plt.tight_layout()
plt.savefig("DCCM_8DD2_plot.png", dpi=300)
print("Analisi completata! Grafico salvato come 'DCCM_8DD2_plot.png'.")