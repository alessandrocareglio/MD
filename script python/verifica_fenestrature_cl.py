import MDAnalysis as mda
import numpy as np

# Carica topologia e traiettoria
u = mda.Universe("step7_production.tpr", "traj_centered.xtc")

# Selezione ioni cloruro e centro di massa proteina
cl_ions = u.select_atoms("name CLA or resname CLA or name CL or resname CL")
protein = u.select_atoms("protein")

print(f"Analisi di {len(cl_ions)} ioni Cl- su {len(u.trajectory)} frame...")

vestibolo_z_min, vestibolo_z_max = 14.0, 22.0
r_vestibolo_max = 12.0

# Registro eventi: {ione_id: {'apicale': int, 'laterale': int}}
eventi_ingresso = {"laterale": 0, "apicale": 0, "dettagli": []}

# Estrazione coordinate centrate sull'asse del recettore
history_pos = []  # buffer temporale per tracciare il vettore di provenienza

for ts in u.trajectory:
    com_xy = protein.center_of_mass()[:2]
    pos = cl_ions.positions.copy()
    pos[:, 0] -= com_xy[0]
    pos[:, 1] -= com_xy[1]
    history_pos.append(pos)

    # Mantieni una finestra di memoria di 15 frame (~400 ps)
    if len(history_pos) > 15:
        history_pos.pop(0)
    else:
        continue

    # Coordinate correnti
    curr_pos = history_pos[-1]
    r_xy = np.sqrt(curr_pos[:, 0] ** 2 + curr_pos[:, 1] ** 2)
    z = curr_pos[:, 2]

    # Identifica ioni attualmente nel vestibolo centrale
    in_vestibule = np.where(
        (z >= vestibolo_z_min)
        & (z <= vestibolo_z_max)
        & (r_xy <= r_vestibolo_max)
    )[0]

    for ion_idx in in_vestibule:
        # Posizione 15 frame prima
        prev_pos = history_pos[0][ion_idx]
        prev_r_xy = np.sqrt(prev_pos[0] ** 2 + prev_pos[1] ** 2)
        prev_z = prev_pos[2]

        # Criterio di provenienza:
        # Se proveniva da R_xy > 25 A a quota simile -> Ingresso Laterale (Fenestratura)
        # Se proveniva da Z > 35 A con R_xy < 15 A -> Ingresso Apicale
        if prev_r_xy > 22.0 and abs(prev_z - z[ion_idx]) < 12.0:
            ion_atom = cl_ions[ion_idx]
            event_id = f"Ion_{ion_atom.id}_Frame_{ts.frame}"
            if event_id not in [e["id"] for e in eventi_ingresso["dettagli"]]:
                eventi_ingresso["laterale"] += 1
                eventi_ingresso["dettagli"].append(
                    {
                        "id": event_id,
                        "tipo": "LATERALE (Fenestratura)",
                        "z": z[ion_idx],
                        "r_prev": prev_r_xy,
                    }
                )

        elif prev_z > 35.0 and prev_r_xy <= 15.0:
            ion_atom = cl_ions[ion_idx]
            event_id = f"Ion_{ion_atom.id}_Frame_{ts.frame}"
            if event_id not in [e["id"] for e in eventi_ingresso["dettagli"]]:
                eventi_ingresso["apicale"] += 1
                eventi_ingresso["dettagli"].append(
                    {
                        "id": event_id,
                        "tipo": "APICALE (Bocca centrale)",
                        "z": z[ion_idx],
                        "r_prev": prev_r_xy,
                    }
                )

print("\n--- RISULTATO TRACCIAMENTO VIE DI ACCESSO CL- ---")
print(f"Ingressi rilevati da FENESTRATURE LATERALI: {eventi_ingresso['laterale']}")
print(f"Ingressi rilevati da BOCCA APICALE:         {eventi_ingresso['apicale']}")

for ev in eventi_ingresso["dettagli"][:10]:
    print(f" -> {ev['id']}: {ev['tipo']} (R_iniziale = {ev['r_prev']:.1f} A)")