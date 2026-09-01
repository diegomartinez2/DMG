#!/usr/bin/env python3
import os
import sys
import gzip
import xml.etree.ElementTree as ET

# Intenta usar la interfaz gráfica para seleccionar archivo; si falla, permite pasar el archivo por argumento.
try:
    import tkinter as tk
    from tkinter import filedialog
    USE_TK = True
except ImportError:
    USE_TK = False

# Mapeo de identificadores de sector a nombres legibles
LOCATION_NAMES = {
    "cluster_500_sector001_macro": "Avarice I",
    "cluster_500_sector002_macro": "Avarice V",
    "cluster_500_sector003_macro": "Avarice IV",
    "cluster_501_macro": "Windfall I",
    "cluster_502_macro": "Windfall III",
    "cluster_503_macro": "Windfall IV",
}

# Coordenadas fijas para ciertas zonas/macros
MACRO_TO_COORDS = {
    "Zone004_Cluster_503_Sector001_macro": "139680,0,-3215.59",
    "Zone003_Cluster_500_Sector003_macro": "-115796.9,0,-96109.38",
    "Zone002_Cluster_503_Sector001_macro": "-2269.199,0,165758.4",
    "Zone006_Cluster_500_Sector002_macro": "-184409.83,0,-3343.0",
    "Zone003_Cluster_501_Sector001_macro": "-97898.44,0,-23632.81",
    "Zone003_Cluster_500_Sector001_macro": "-97898.44,0,-23632.81",
    "Zone005_Cluster_501_Sector001_macro": "-23023.44,0,-160132.8",
    "Zone004_Cluster_504_Sector001_macro": "157616.8,0,-84263.63",
    "tzoneCluster_500_Sector002SHCon9_GateZone_macro": "-39827.2,0,200000",
    "tzoneCluster_500_Sector002SHCon5_GateZone_macro": "-50000,0,200000",
    "Zone001_Cluster_500_Sector003_macro": "-117046.9,0,82593.75",
    "Zone001_Cluster_501_Sector001_macro": "-127062.5,0,133718.8",
    "Zone005_Cluster_504_Sector001_macro": "-23023.44,0,-160132.8",
    "Zone004_Cluster_500_Sector001_macro": "103164.1,0,-32492.19",
    "Zone002_Cluster_500_Sector001_macro": "-10960.94,0,112070.3",
    "Zone003_Cluster_500_Sector002_macro": "-7132.69,0,531.25",
    "Zone002_Cluster_500_Sector002_macro": "-118489.4,0,-131687.5",
    "tzoneCluster_500_Sector003SHCon9_GateZone_macro": "50000,0,200000",
    "Zone005_Cluster_500_Sector002_macro": "98694.83,0,-126406.3",
    "tzoneCluster_500_Sector001SHCon2_GateZone_macro": "66960.9,0,-207957",
    "tzoneCluster_500_Sector003SHCon5_GateZone_macro": "40054.7,0,200000",
    "tzoneCluster_500_Sector001SHCon10_GateZone_macro": "-45738.28,0,-200000",
    "Zone005_Cluster_500_Sector003_macro": "34671.88,0,-152515.6",
    "Zone001_Cluster_504_Sector001_macro": "-127062.5,0,133718.8",
    "Zone003_Cluster_502_Sector001_macro": "-97898.44,0,-23632.81",
    "Zone002_Cluster_504_Sector001_macro": "137904.1,0,135199.3",
    "Zone005_Cluster_503_Sector001_macro": "-23023.44,0,-160132.8",
    "Zone001_Cluster_500_Sector001_macro": "-127062.5,0,133718.8",
    "Zone004_Cluster_501_Sector001_macro": "103164.1,0,-32492.19",
    "tzoneCluster_500_Sector001SHCon4_GateZone_macro": "76886.7,0,-207957",
    "Zone006_Cluster_502_Sector001_macro": "47000.0,0,172000.0",
    "Zone002_Cluster_500_Sector003_macro": "44015.63,0,158125",
    "tzoneCluster_500_Sector001SHCon6_GateZone_macro": "-55785.2,0,-200000",
    "Zone002_Cluster_501_Sector001_macro": "-10960.94,0,112070.3",
    "Zone001_Cluster_500_Sector002_macro": "-78593.63,0,157937.5",
    "Zone003_Cluster_503_Sector001_macro": "-199365.3,0,-55454.45",
    "Zone001_Cluster_502_Sector001_macro": "-127062.5,0,133718.8",
    "Zone003_Cluster_504_Sector001_macro": "-226443.3,0,-97194.67",
    "Zone005_Cluster_500_Sector001_macro": "-23023.44,0,-160132.8",
    "Zone004_Cluster_500_Sector003_macro": "110359.4,0,7250",
    "Zone001_Cluster_503_Sector001_macro": "-190039.3,0,180073.9",
    "Zone005_Cluster_502_Sector001_macro": "-23023.44,0,-160132.8",
    "Zone004_Cluster_500_Sector002_macro": "101721.7,0,100656.3",
    "Zone002_Cluster_502_Sector001_macro": "108009.3,0,115933.7",
    "Zone004_Cluster_502_Sector001_macro": "103164.1,0,-32492.19"
}

def seleccionar_archivo():
    if len(sys.argv) > 1:
        return sys.argv[1]

    if USE_TK:
        root = tk.Tk()
        root.withdraw()
        path_default = os.path.expanduser("~/.local/share/Egosoft/X4")
        filename = filedialog.askopenfilename(
            initialdir=path_default if os.path.exists(path_default) else os.path.expanduser("~"),
            title="Selecciona la partida de X4",
            filetypes=[("Partida X4", "*.xml *.xml.gz"), ("Todos", "*.*")]
        )
        return filename
    else:
        return input("Introduce la ruta al archivo de partida (.xml / .xml.gz): ").strip()

def obtener_nodos_con_padres(tree):
    parent_map = {c: p for p in tree.iter() for c in p}
    return parent_map

def calcular_posicion(node, parent_map):
    total_x, total_y, total_z = 0.0, 0.0, 0.0
    ubicacion = None

    curr = node
    while curr is not None:
        cls = curr.attrib.get('class', '')
        macro = curr.attrib.get('macro', '')

        if cls == "galaxy":
            break

        if cls == "sector" and not ubicacion:
            ubicacion = LOCATION_NAMES.get(macro, macro)
        elif cls == "cluster" and not ubicacion:
            ubicacion = LOCATION_NAMES.get(macro, macro)

        x, y, z = 0.0, 0.0, 0.0

        if macro in MACRO_TO_COORDS:
            parts = [float(val) for val in MACRO_TO_COORDS[macro].split(',')]
            x, y, z = parts[0], parts[1], parts[2]
        else:
            offset = curr.find('offset')
            if offset is not None:
                pos = offset.find('position')
                if pos is not None:
                    x = float(pos.attrib.get('x', 0))
                    y = float(pos.attrib.get('y', 0))
                    z = float(pos.attrib.get('z', 0))

        total_x += x
        total_y += y
        total_z += z

        curr = parent_map.get(curr)

    return ubicacion, total_x / 1000.0, total_y / 1000.0, total_z / 1000.0

def main():
    filepath = seleccionar_archivo()
    if not filepath or not os.path.exists(filepath):
        print("No se seleccionó ningún archivo válido.")
        return

    print(f"Cargando archivo: {filepath}")
    print("Analizando XML (esto puede usar bastante memoria RAM)...")

    if filepath.endswith('.gz'):
        with gzip.open(filepath, 'rb') as f:
            tree = ET.parse(f)
    else:
        tree = ET.parse(filepath)

    parent_map = obtener_nodos_con_padres(tree)
    root = tree.getroot()

    print("\nBuscando Bóvedas de Datos de Erlking...\n")
    resultados = []

    for comp in root.findall(".//component"):
        macro = comp.attrib.get('macro', '')
        if "landmarks_erlking_vault" in macro:
            ubicacion, x_km, y_km, z_km = calcular_posicion(comp, parent_map)
            code = comp.attrib.get('code', 'N/A')

            # Buscar planos dentro del nodo
            blueprints = []
            bp_nodes = comp.findall(".//blueprints")
            for bp in bp_nodes:
                if bp.text:
                    blueprints.append(bp.text.strip())

            resultados.append({
                "Ubicación": ubicacion,
                "X (km)": round(x_km, 2),
                "Y (km)": round(y_km, 2),
                "Z (km)": round(z_km, 2),
                "Código": code,
                "Bóveda": macro
            })

    # Mostrar resultados en consola
    if not resultados:
        print("No se encontraron bóvedas del Erlking en esta partida.")
    else:
        header = f"{'Ubicación':<18} | {'X (km)':<10} | {'Y (km)':<10} | {'Z (km)':<10} | {'Código':<10} | {'Bóveda'}"
        print("-" * len(header))
        print(header)
        print("-" * len(header))
        for r in resultados:
            print(f"{r['Ubicación']:<18} | {r['X (km)']:<10} | {r['Y (km)']:<10} | {r['Z (km)']:<10} | {r['Código']:<10} | {r['Bóveda']}")
        print("-" * len(header))

if __name__ == "__main__":
    main()
