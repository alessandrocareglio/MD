import MDAnalysis as mda
import numpy as np

traiettoria = "gaba_30ns_noPBC.xtc"
topologia = "struttura_30ns.pdb"

u = mda.Universe(topologia, traiettoria)
n_frames = len(u.trajectory)
n_atoms = u.atoms.n_atoms
print(f"Traiettoria caricata: {n_atoms} atomi, {n_frames} frame.")

# Identificazione per indici atomici:
# Gli ultimi 73 atomi sono i ligandi (43 ATP + 30 GABA)
# In base all'ordinamento standard di GROMACS/CHARMM-GUI:
# Proteina/Core: atomi da 0 a 34028 (34029 atomi)
# ATP: 43 atomi (da 34029 a 34071)
# GABA: 30 atomi (da 34072 a 34101)

idx_core = np.arange(0, 34029)
idx_atp = np.arange(34029, 34072)

core_atoms = u.atoms[idx_core]
atp_atoms = u.atoms[idx_atp]

print(f"Atomi core proteico per allineamento: {len(core_atoms)}")
print(f"Atomi ATP selezionati: {len(atp_atoms)}")

# Riferimento: frame 0
u.trajectory[0]
ref_core_pos = core_atoms.positions.copy()
ref_core_com = ref_core_pos.mean(axis=0)
ref_core_centered = ref_core_pos - ref_core_com

ref_atp_pos = atp_atoms.positions.copy()
ref_atp_centered = ref_atp_pos - ref_core_com

rmsd_list = []

print("Calcolo RMSD dell'ATP in corso...")
for ts in u.trajectory:
    curr_core_pos = core_atoms.positions
    curr_core_com = curr_core_pos.mean(axis=0)

    # Allineamento Kabsch sul core proteico
    P = curr_core_pos - curr_core_com
    Q = ref_core_centered
    H = np.dot(P.T, Q)
    V, S, Wt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(np.dot(Wt.T, V.T)))
    U_mat = np.dot(Wt.T, np.array([[1, 0, 0], [0, 1, 0], [0, 0, d]]))
    U_rot = np.dot(U_mat, V.T)

    # Rototraslazione delle coordinate dell'ATP
    curr_atp_centered = atp_atoms.positions - curr_core_com
    curr_atp_rot = np.dot(curr_atp_centered, U_rot.T)

    # RMSD rispetto alla posa di partenza (frame 0)
    diff = curr_atp_rot - ref_atp_centered
    rmsd = np.sqrt(np.mean(np.sum(diff**2, axis=1)))
    rmsd_list.append([ts.frame, ts.time, rmsd])

rmsd_arr = np.array(rmsd_list)

with open("rmsd_atp_30ns.dat", "w") as f:
    f.write("# Frame  Tempo(ps)  RMSD_ATP(Angstrom)\n")
    for r in rmsd_arr:
        f.write(f"{int(r[0]):6d}  {r[1]:10.2f}  {r[2]:10.4f}\n")

print("\n--- RISULTATO DINAMICA ATP (30 ns in soluzione) ---")
print(
    f"RMSD medio ATP: {np.mean(rmsd_arr[:, 2]):.3f} +/-"
    f" {np.std(rmsd_arr[:, 2]):.3f} Angstrom"
)
print(f"RMSD max:       {np.max(rmsd_arr[:, 2]):.3f} Angstrom")
print(f"RMSD min:       {np.min(rmsd_arr[:, 2]):.3f} Angstrom")
print("Dati esportati con successo in: rmsd_atp_30ns.dat")