import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os # Importar para manejar rutas de archivos

# --- Parámetros de la medición (DEBES MODIFICAR ESTOS VALORES) ---
# Diámetro del bucle transmisor (en metros)
DIAMETRO_LOOP = 40.0
# Corriente del bucle transmisor (en Amperios)
CORRIENTE_TX = 5.0
# Número de vueltas en la bobina receptora
NUM_VUELTAS_RX = 1

# --- Constantes físicas ---
MU_0 = 4 * np.pi * 1e-7  # Permeabilidad magnética del vacío (H/m)

# Calcular el área efectiva de la bobina receptora (Asume que es igual al Tx, por simplicidad)
AREA_RX = np.pi * (DIAMETRO_LOOP / 2)**2 * NUM_VUELTAS_RX
# Momento dipolar del transmisor
MOMENTO_DIPOLAR = CORRIENTE_TX * (np.pi * (DIAMETRO_LOOP / 2)**2)
# Constante para el cálculo de la Resistividad Aparente (Aproximación Late-Time)
C_TEM = 1 / (20 * np.pi) * (4 * np.pi * MU_0)**(5 / 2)


def analizar_datos_tem(ruta_archivo):
    """
    Carga datos de TEM, calcula la resistividad aparente y genera una gráfica.
    """
    try:
        # 1. Cargar datos
        df = pd.read_csv(ruta_archivo)

        # Verificar columnas
        if 'Tiempo (s)' not in df.columns or 'Voltaje (V)' not in df.columns:
            print("Error: El archivo CSV debe contener las columnas 'Tiempo (s)' y 'Voltaje (V)'.")
            return

        tiempo = df['Tiempo (s)'].values
        voltaje = df['Voltaje (V)'].values

        # 2. Calcular la Resistividad Aparente (ρa)

        # Voltaje normalizado por el área efectiva del receptor
        voltaje_normalizado = np.abs(voltaje) / AREA_RX

        # Aplicar la fórmula de Resistividad Aparente (Aproximación Late-Time):
        # rho_a = ( (C_TEM * M) / V_t_normalizado )^(2/3) * t^(5/3)
        # Se usa np.abs() para asegurar que la división sea positiva, ya que logarítmicamente
        # y físicamente el voltaje es siempre positivo en el modelo de difusión.
        rho_aparente = (C_TEM * MOMENTO_DIPOLAR / voltaje_normalizado)**(2/3) * (tiempo)**(5/3)

        df['Resistividad Aparente (Ohm.m)'] = rho_aparente

        # 3. Graficar los resultados

        fig, ax1 = plt.subplots(figsize=(10, 7))

        # Eje Y izquierda: Voltaje (Señal de Decaimiento)
        color = 'tab:red'
        ax1.set_xlabel('Tiempo (s) $\\rightarrow$ Profundidad', fontsize=12)
        ax1.set_ylabel('Voltaje Inducido $|\mathbf{V}|$(V)', color=color, fontsize=12)
        ax1.loglog(tiempo, np.abs(voltaje), marker='.', linestyle='-', color=color, label='Voltaje Medido')
        ax1.tick_params(axis='y', labelcolor=color)
        ax1.grid(True, which="both", ls="--", alpha=0.5)

        # Configurar escalas logarítmicas e invertir el eje X
        ax1.set_xscale('log')
        ax1.set_yscale('log')
        # Invertir el eje X para que los "tiempos tempranos" (superficie) queden a la derecha (opcional, pero común)
        # ax1.invert_xaxis()

        # Eje Y derecha: Resistividad Aparente
        ax2 = ax1.twinx()
        color = 'tab:blue'
        ax2.set_ylabel('Resistividad Aparente $\\rho_a$ ($\Omega \cdot$m)', color=color, fontsize=12)
        ax2.loglog(tiempo, rho_aparente, marker='o', linestyle='none', color=color, label='Resistividad Aparente')
        ax2.tick_params(axis='y', labelcolor=color)
        ax2.set_yscale('log')

        plt.title(f'Curva de Decaimiento TEM y Resistividad Aparente (Tx: {DIAMETRO_LOOP}m, {CORRIENTE_TX}A)', fontsize=14)
        fig.tight_layout()
        plt.show()

        print("\n--- Resultados del Análisis Preliminar (Primeras 5 filas) ---")
        print(df[['Tiempo (s)', 'Voltaje (V)', 'Resistividad Aparente (Ohm.m)']].head())
        print(f"\nResistividad Aparente Mínima: {rho_aparente.min():.2f} Ohm.m")
        print(f"Resistividad Aparente Máxima: {rho_aparente.max():.2f} Ohm.m")

    except FileNotFoundError:
        print(f"Error: No se encontró el archivo en la ruta especificada: {ruta_archivo}")
    except Exception as e:
        print(f"Ocurrió un error durante el procesamiento: {e}. Verifique sus parámetros de medición.")

# --- PUNTO DE ENTRADA PRINCIPAL (MAIN) ---
if __name__ == "__main__":
    # 1. Crear un archivo CSV de prueba para que el script sea ejecutable
    # Estos son datos típicos de decaimiento (el voltaje cae con el tiempo)
    datos_prueba = {
        'Tiempo (s)': np.logspace(-4, -1, 50), # Tiempos desde 0.1ms hasta 100ms
        'Voltaje (V)': 1e-3 * np.exp(-1000 * np.logspace(-4, -1, 50)) + 1e-6 # Señal con decaimiento
    }
    df_prueba = pd.DataFrame(datos_prueba)

    # Crear un archivo temporal de prueba en la misma carpeta
    nombre_archivo_prueba = 'datos_tem_ejemplo.csv'
    df_prueba.to_csv(nombre_archivo_prueba, index=False)
    print(f"Archivo de prueba '{nombre_archivo_prueba}' creado exitosamente.")

    # 2. Ejecutar la función con el archivo de prueba
    print("\nEjecutando análisis TEM con datos de ejemplo...")
    analizar_datos_tem(nombre_archivo_prueba)

    # 3. Limpiar el archivo de prueba (opcional)
    # os.remove(nombre_archivo_prueba)
    # print(f"\nArchivo de prueba '{nombre_archivo_prueba}' eliminado.")
