import MDAnalysis as mda
import numpy as np

traiettoria = "gaba_30ns_noPBC.xtc"
topologia = "struttura_30ns.pdb"

u = mda.Universe(topologia, traiettoria)
n_frames = len(u.trajectory)
print(f"Caricati {u.atoms.n_atoms} atomi su {n_frames} frame.")

# 1. VERIFICA COMPOSIZIONE
waters = u.select_atoms("name OH2 or resname TIP3 or resname SOL")
cl_ions = u.select_atoms("name CLA or resname CLA or name CL or resname CL")
print(f"Acque rilevate: {len(waters)}")
print(f"Ioni Cloro rilevati: {len(cl_ions)}")

# 2. SE PRESENTI ACQUE, CALCOLA IL PROFILO
if len(waters) > 0:
    sel_center = u.select_atoms("protein")
    cylinder_radius = 5.0
    z_bins = np.arange(-40, 41, 1.0)
    water_counts = np.zeros(len(z_bins) - 1)
    frames_count = 0

    for ts in u.trajectory[::5]:
        frames_count += 1
        com_xy = sel_center.center_of_mass()[:2]
        pos = waters.positions
        dx = pos[:, 0] - com_xy[0]
        dy = pos[:, 1] - com_xy[1]
        r = np.sqrt(dx**2 + dy**2)
        inside = pos[r <= cylinder_radius]
        counts, _ = np.histogram(inside[:, 2], bins=z_bins)
        water_counts += counts

    dz = 1.0
    r_nm = cylinder_radius / 10.0
    vol_nm3 = np.pi * (r_nm**2) * (dz / 10.0)
    density = (water_counts / frames_count) / vol_nm3
    z_centers = 0.5 * (z_bins[:-1] + z_bins[1:])

    with open("densita_acqua_poro_30ns.dat", "w") as f:
        f.write("# Quota_Z(A)  Densita_30ns(molecole/nm3)\n")
        for z, d in zip(z_centers, density):
            f.write(f"{z:10.2f}  {d:12.4f}\n")
    print("Densità acqua scritta in: densita_acqua_poro_30ns.dat")

# 3. SE PRESENTI CLORURI, SALVA LE QUOTE Z
if len(cl_ions) > 0:
    z_traj = np.zeros((n_frames, len(cl_ions)))
    for idx, ts in enumerate(u.trajectory):
        z_traj[idx, :] = cl_ions.positions[:, 2]
    np.save("z_trajectories_cl_30ns.npy", z_traj)
    print("Traiettorie Cloro salvate in: z_trajectories_cl_30ns.npy")