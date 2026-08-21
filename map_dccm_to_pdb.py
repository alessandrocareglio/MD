import MDAnalysis as mda
from MDAnalysis.analysis import align
import numpy as np

# 1. Carica universo
u = mda.Universe("prod.tpr", "prod.xtc")

# Allinea la traiettoria sui C-alpha
calphas = u.select_atoms("protein and name CA")
ref = u.select_atoms("protein and name CA")
align.AlignTraj(u, ref, select="protein and name CA", in_memory=True).run()

n_ca = len(calphas)
print(f"Atomi C-alfa analizzati: {n_ca}")

# 2. Identifica i C-alpha vicini all'ATP (< 7 Å)
nearby_ca_atoms = u.select_atoms("protein and name CA and around 7.0 resname ATP")
ca_indices_in_universe = list(calphas.indices)
atp_site_indices = [ca_indices_in_universe.index(a.index) for a in nearby_ca_atoms]
print(f"Trovati {len(atp_site_indices)} residui C-alfa nel sito di legame dell'ATP.")

# 3. Estrai coordinate e calcola matrice DCCM
coords = []
for ts in u.trajectory:
    coords.append(calphas.positions)
coords = np.array(coords)

mean_pos = np.mean(coords, axis=0)
dr = coords - mean_pos

dot_prod = np.zeros((n_ca, n_ca))
for frame in range(dr.shape[0]):
    dot_prod += np.dot(dr[frame], dr[frame].T)

mean_dot = dot_prod / dr.shape[0]
diag = np.diag(mean_dot)
dccm = mean_dot / np.sqrt(np.outer(diag, diag))

# 4. Calcola la correlazione media con il sito ATP
atp_corr = np.mean(dccm[:, atp_site_indices], axis=1)

# Mappa per resindex (univoco per ogni residuo nella topologia)
corr_dict = {calphas[i].resindex: atp_corr[i] * 100.0 for i in range(n_ca)}

# 5. Inizializza l'attributo tempfactors nell'universo per tutti gli atomi
protein = u.select_atoms("protein")
all_tempfactors = np.zeros(u.atoms.n_atoms)

for atom in protein:
    all_tempfactors[atom.index] = corr_dict.get(atom.resindex, 0.0)

# Aggiungi formalmente l'attributo B-factor alla topologia
u.add_TopologyAttr('tempfactors', all_tempfactors)

# 6. Salva la proteina con i B-factors mappati
output_pdb = "8DD2_pathway_correlated.pdb"
protein.write(output_pdb)
print(f"File salvato con successo: '{output_pdb}'")