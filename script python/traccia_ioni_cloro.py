import MDAnalysis as mda
import numpy as np

# Carica topologia e traiettoria centrata
u = mda.Universe("step7_production.tpr", "traj_centered.xtc")

# Seleziona tutti gli ioni Cloro
cl_ions = u.select_atoms("name CLA or resname CLA or name CL or resname CL")
n_ions = len(cl_ions)
print(f"Numero totale di ioni Cloro nel sistema: {n_ions}")

# Array per memorizzare le coordinate Z di tutti gli ioni nel tempo
n_frames = len(u.trajectory)
z_trajectories = np.zeros((n_frames, n_ions))

print("Scansione della traiettoria in corso...")
for frame_idx, ts in enumerate(u.trajectory):
    z_trajectories[frame_idx, :] = cl_ions.positions[:, 2]

# Analisi delle traslocazioni e attraversamenti
z_min = np.min(z_trajectories, axis=0)
z_max = np.max(z_trajectories, axis=0)
delta_z = z_max - z_min

# Identifica ioni che hanno oltrepassato il centro della membrana (Z = 0)
penetrated_core = np.where(z_min < 0.0)[0]

# Identifica ioni che hanno attraversato l'intero bilayer (es. da Z > +20 a Z < -20)
full_crossing = np.where((z_max > 20.0) & (z_min < -20.0))[0]

print("\n--- RISULTATI DEL CONTROLLO IONI ---")
print(f"Ioni che sono scesi sotto Z = 0 A: {len(penetrated_core)}")
if len(penetrated_core) > 0:
    for idx in penetrated_core:
        atom = cl_ions[idx]
        print(f"  -> Ione ID {atom.id} (Residuo {atom.resid}): Z_min = {z_min[idx]:.2f} A, Z_max = {z_max[idx]:.2f} A")

print(f"\nIoni con traslocazione trans-membrana completa (Z > +20 -> Z < -20): {len(full_crossing)}")
if len(full_crossing) > 0:
    for idx in full_crossing:
        atom = cl_ions[idx]
        print(f"  -> TRASLOCAZIONE RILEVATA: Ione ID {atom.id} (Residuo {atom.resid})")

# Salva le traiettorie Z per generare un grafico se necessario
np.save("z_trajectories_cl.npy", z_trajectories)
print("\nDati completi salvati in z_trajectories_cl.npy")