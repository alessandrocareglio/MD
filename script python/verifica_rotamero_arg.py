import MDAnalysis as mda
import numpy as np

# Carica topologia e traiettoria
u = mda.Universe("step7_production.tpr", "traj_centered.xtc")

# Selezione specifica dell'arginina identificata
res_target = u.select_atoms("segid seg_0_PROA and resid 263").residues[0]
print(f"Residuo selezionato: {res_target.resname} {res_target.resid} ({res_target.segid})")

# Selezione atomi per Chi1 (N - CA - CB - CG)
N  = res_target.atoms.select_atoms("name N")[0]
CA = res_target.atoms.select_atoms("name CA")[0]
CB = res_target.atoms.select_atoms("name CB")[0]
CG = res_target.atoms.select_atoms("name CG")[0]
CD = res_target.atoms.select_atoms("name CD")[0]

gruppo_chi1 = mda.AtomGroup([N, CA, CB, CG])
gruppo_chi2 = mda.AtomGroup([CA, CB, CG, CD])

tempi = []
chi1_val = []
chi2_val = []

print("Calcolo torsioni Chi1 e Chi2 sui 5001 frame in corso...")
for ts in u.trajectory[::5]:  # Campiona ogni 5 frame (1000 punti lungo i 100 ns)
    tempi.append(ts.time / 1000.0)  # Convertito in ns
    chi1_val.append(gruppo_chi1.dihedral.value())
    chi2_val.append(gruppo_chi2.dihedral.value())

tempi = np.array(tempi)
chi1_val = np.array(chi1_val)
chi2_val = np.array(chi2_val)

# Salva i dati numerici
with open("dihedral_arg263.dat", "w") as f:
    f.write("# Tempo(ns)   Chi1(deg)   Chi2(deg)\n")
    for t, c1, c2 in zip(tempi, chi1_val, chi2_val):
        f.write(f"{t:10.3f}   {c1:10.2f}   {c2:10.2f}\n")

print("\n--- RISULTATO DINAMICA CATENA LATERALE (ARG 263) ---")
print(f"Chi1 medio: {np.mean(chi1_val):.1f}° +/- {np.std(chi1_val):.1f}° (Range: {np.min(chi1_val):.1f}° -> {np.max(chi1_val):.1f}°, Escursione: {np.ptp(chi1_val):.1f}°)")
print(f"Chi2 medio: {np.mean(chi2_val):.1f}° +/- {np.std(chi2_val):.1f}° (Range: {np.min(chi2_val):.1f}° -> {np.max(chi2_val):.1f}°, Escursione: {np.ptp(chi2_val):.1f}°)")

if np.ptp(chi1_val) > 60 or np.ptp(chi2_val) > 60:
    print("\n-> TRANSIZIONE CONFERMATA: La catena laterale di ARG 263 ha compiuto rotazioni rotameriche reali durante la dinamica.")
else:
    print("\n-> ROTAZIONE NON OSSERVATA: Il residuo è rimasto fermo nello stesso pozzo conformazionale.")