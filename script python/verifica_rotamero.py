import MDAnalysis as mda
import numpy as np

# Topologia e traiettoria di produzione
u = mda.Universe("step7_production.tpr", "traj_centered.xtc")

# Inserisci il numero del residuo osservato nel morph
res_num = 156  # Sostituisci con il numero reale del residuo (es. 156 o altro)
residuo = u.select_atoms(f"resid {res_num}")
print(f"Residuo analizzato: {residuo.resnames[0]} {res_num}")

# Atomi per il diedro chi1 (N - CA - CB - CG)
try:
    atom_N = u.select_atoms(f"resid {res_num} and name N")[0]
    atom_CA = u.select_atoms(f"resid {res_num} and name CA")[0]
    atom_CB = u.select_atoms(f"resid {res_num} and name CB")[0]
    atom_CG = u.select_atoms(f"resid {res_num} and (name CG or name CG1)")[0]
except IndexError:
    print("Impossibile trovare la quaterna N-CA-CB-CG per questo residuo.")
    exit(1)

diedro_group = mda.core.groups.AtomGroup([atom_N, atom_CA, atom_CB, atom_CG])

tempi = []
angoli = []

print("Calcolo torsione chi1 sui 5001 frame in corso...")
for ts in u.trajectory[::5]:  # campiona ogni 5 frame (1000 punti)
    # Calcolo diedro in gradi
    angle = diedro_group.dihedral.value()
    tempi.append(ts.time / 1000.0)  # in ns
    angoli.append(angle)

tempi = np.array(tempi)
angoli = np.array(angoli)

with open("dihedral_residuo.dat", "w") as f:
    f.write("# Tempo(ns)   Chi1_Angle(deg)\n")
    for t, a in zip(tempi, angoli):
        f.write(f"{t:10.3f}   {a:10.2f}\n")

print("\n--- RISULTATO DINAMICA CATENA LATERALE ---")
print(f"Angolo medio: {np.mean(angoli):.1f}° +/- {np.std(angoli):.1f}°")
print(f"Escursione angolare: da {np.min(angoli):.1f}° a {np.max(angoli):.1f}°")
if np.ptp(angoli) > 60:
    print("-> TRANSIZIONE CONFERMATA: La catena laterale ha ruotato realmente nella traiettoria MD.")
else:
    print("-> ROTAZIONE NON OSSERVATA: Il residuo è rimasto fermo nello stesso rotamero (il morph era un artefatto di interpolazione).")