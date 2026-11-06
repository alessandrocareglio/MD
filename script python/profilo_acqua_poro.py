import MDAnalysis as mda
import numpy as np

# Carica topologia e traiettoria centrata
u = mda.Universe("step7_production.tpr", "traj_centered.xtc")

# Trova il centro di massa del canale M2 per centrare il cilindro
# Nel GABA_A M2 corrisponde approssimativamente ai residui attorno al gate
# (usiamo l'asse del canale prendendo la proteina transmembrana)
sel_tmd = u.select_atoms("protein and (resid 250:300)")

# Raggio del cilindro attorno al poro in Angstrom (5 A racchiude il lume interno)
cylinder_radius = 5.0

# Bin lungo Z da -40 A a +40 A (passo di 1 A)
z_bins = np.arange(-40, 41, 1.0)
water_counts = np.zeros(len(z_bins) - 1)
n_frames = 0

print("Calcolo densità dell'acqua nel cilindro del poro in corso...")

for ts in u.trajectory[::10]:  # campiona ogni 10 frame per rapidità
    n_frames += 1
    # Coordinate del centro del poro sul piano XY
    com_xy = sel_tmd.center_of_mass()[:2]

    # Seleziona gli ossigeni dell'acqua
    waters = u.select_atoms("resname TIP3 and name OH2")
    pos = waters.positions

    # Distanza radiale dall'asse del canale
    dx = pos[:, 0] - com_xy[0]
    dy = pos[:, 1] - com_xy[1]
    r = np.sqrt(dx**2 + dy**2)

    # Filtra solo le acque dentro il cilindro del poro
    inside_pore = pos[r <= cylinder_radius]

    # Istogramma lungo l'asse Z
    counts, _ = np.histogram(inside_pore[:, 2], bins=z_bins)
    water_counts += counts

# Calcolo del volume di ciascun segmento cilindrico (in nm^3)
dz = 1.0  # Angstrom = 0.1 nm
r_nm = cylinder_radius / 10.0  # 0.5 nm
volume_bin_nm3 = np.pi * (r_nm**2) * (dz / 10.0)

# Densità media (molecole / nm^3)
water_density = (water_counts / n_frames) / volume_bin_nm3
z_centers = 0.5 * (z_bins[:-1] + z_bins[1:])

# Salva i dati su file
with open("densita_acqua_poro_cilindro.dat", "w") as f:
    f.write("# Quota_Z(A)  Densita_Acqua_Poro(molecole/nm3)\n")
    for z, dens in zip(z_centers, water_density):
        f.write(f"{z:10.2f}  {dens:12.4f}\n")

print(
    "Calcolo completato! Risultati scritti in densita_acqua_poro_cilindro.dat"
)