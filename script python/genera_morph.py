import MDAnalysis as mda
import numpy as np

# 1. Carica topologia e traiettoria
topologia = "step7_production.tpr"
traiettoria = "traj_centered.xtc"

u = mda.Universe(topologia, traiettoria)

# 2. Definisci i frame di riferimento (indicativi: frame 0 per chiuso, frame 940 per picco a ~18.8 ns)
# (Ciascun frame corrisponde a 20 ps: 18800 ps / 20 ps = frame 940)
frame_chiuso = 0
frame_aperto = 940

# Isola solo la proteina (o la proteina + ligandi) per mantenere il file compatto e leggero
selezione = u.select_atoms("protein")

# Estrai coordinate stato chiuso
u.trajectory[frame_chiuso]
pos_chiuso = selezione.positions.copy()

# Estrai coordinate stato aperto
u.trajectory[frame_aperto]
pos_aperto = selezione.positions.copy()

# 3. Parametri del morphing
n_passi = 40  # Numero di frame intermedi generati
alpha_valori = np.linspace(0.0, 1.0, n_passi)

print(
    f"Generazione morphing lineare tra frame {frame_chiuso} e frame"
    f" {frame_aperto}..."
)

# 4. Scrittura del file multi-frame PDB
nome_output = "morph_poro_lineare.pdb"
with mda.Writer(nome_output, selezione.n_atoms) as W:
    for alpha in alpha_valori:
        # Interpolazione lineare pura
        pos_interp = (1.0 - alpha) * pos_chiuso + alpha * pos_aperto
        selezione.positions = pos_interp
        W.write(selezione)

print(f"File generato con successo: {nome_output} ({n_passi} frame)")