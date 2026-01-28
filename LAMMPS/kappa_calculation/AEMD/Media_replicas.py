#!/usr/local/bin/python
import os
import pandas as pd
import numpy as np
import glob

def process_results():
    base_dir = "RESULTS"
    # Buscamos todas las combinaciones de parámetros (exceptuando la carpeta REP)
    # Estructura: RESULTS/SIZE_*/DT_*/TS_*/
    case_patterns = glob.glob(os.path.join(base_dir, "SIZE_*", "DT_*", "TS_*"))

    for case_path in case_patterns:
        all_data = []

        # Buscar todas las réplicas dentro de este caso
        replica_dirs = glob.glob(os.path.join(case_path, "REP_*"))

        if not replica_dirs:
            continue

        for rep_dir in replica_dirs:
            file_path = os.path.join(rep_dir, "delta_T_raw.dat")
            if os.path.exists(file_path):
                # LAMMPS ave/time tiene 2 líneas de comentario iniciales
                # El formato suele ser: TimeStep Number-of-bins v_delta_T
                try:
                    df = pd.read_csv(file_path, sep='\s+', skiprows=2, names=['Step', 'Count', 'DeltaT'])
                    all_data.append(df.set_index('Step')['DeltaT'])
                except Exception as e:
                    print(f"Error leyendo {file_path}: {e}")

        if all_data:
            # Concatenar réplicas por el índice (Step)
            combined = pd.concat(all_data, axis=1)

            # Calcular media y desviación estándar
            mean_series = combined.mean(axis=1)
            std_series = combined.std(axis=1)

            # Crear DataFrame de salida
            final_df = pd.DataFrame({
                'Mean_DeltaT': mean_series,
                'Std_DeltaT': std_series,
                'N_Replicas': combined.count(axis=1)
            }).reset_index()

            # Guardar en la carpeta del caso (fuera de REP_X)
            output_name = os.path.join(case_path, "averaged_delta_T.dat")
            final_df.to_csv(output_name, sep='\t', index=False)
            print(f"Generado: {output_name}")

if __name__ == "__main__":
    process_results()
