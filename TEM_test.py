import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

'''
script en Python que realiza dos pasos esenciales para el análisis preliminar:

    Lectura y Visualización: Carga los datos de Tiempo vs. Voltaje (la señal de decaimiento).

    Cálculo de la Resistividad Aparente: Convierte el voltaje en una métrica geofísica inicial (ρa​) y lo grafica en una escala logarítmica.

Para que el script funcione, debes asegurarte de que tu archivo CSV contenga las columnas Tiempo (s) y Voltaje (V)
'''

# --- Parámetros de la medición (DEBES MODIFICAR ESTOS VALORES) ---
# Diámetro del bucle transmisor (en metros)
DIAMETRO_LOOP = 40.0
# Corriente del bucle transmisor (en Amperios)
CORRIENTE_TX = 5.0
# Número de vueltas en la bobina receptora
NUM_VUELTAS_RX = 1

# --- Constantes físicas ---
MU_0 = 4 * np.pi * 1e-7  # Permeabilidad magnética del vacío (H/m)
AREA_RX = np.pi * (DIAMETRO_LOOP / 2)**2 * NUM_VUELTAS_RX  # Área efectiva del receptor

def analizar_datos_tem(ruta_archivo):
    """
    Carga datos de TEM, calcula la resistividad aparente y genera una gráfica de decaimiento.

    Args:
        ruta_archivo (str): Ruta al archivo CSV con columnas 'Tiempo (s)' y 'Voltaje (V)'.
    """
    try:
        # 1. Cargar datos
        df = pd.read_csv(ruta_archivo)

        # Verificar que las columnas existan
        if 'Tiempo (s)' not in df.columns or 'Voltaje (V)' not in df.columns:
            print("Error: El archivo CSV debe contener las columnas 'Tiempo (s)' y 'Voltaje (V)'.")
            return

        tiempo = df['Tiempo (s)'].values
        voltaje = df['Voltaje (V)'].values

        # 2. Calcular la Resistividad Aparente (ρa)

        # Fórmula simplificada para el cálculo de la resistividad aparente (Late-Time Approximation)
        # Esta aproximación es más precisa para los tiempos tardíos (mayor profundidad).

        # Ecuación de la aproximación de 'Late-Time' (resistividad aparente)
        # rho_a = ( (C * M) / V_t )^(2/3) * t^(5/3)
        # Donde:
        # M = Momento dipolar (Corriente * Area_TX)
        # V_t = Voltaje de decaimiento
        # C = Constante que contiene 1/(20*pi) * (4*pi*mu_0)^(5/2)

        MOMENTO_DIPOLAR = CORRIENTE_TX * (np.pi * (DIAMETRO_LOOP / 2)**2)
        C = 1 / (20 * np.pi) * (4 * np.pi * MU_0)**(5 / 2)

        # Voltaje normalizado por el área del receptor
        voltaje_normalizado = voltaje / AREA_RX

        # Aplicar la fórmula de resistividad aparente (rho_a)
        # Se usa np.abs para manejar voltajes negativos (aunque se prefiere evitar en este cálculo)
        rho_aparente = (C * MOMENTO_DIPOLAR / np.abs(voltaje_normalizado))**(2/3) * (tiempo)**(5/3)

        df['Resistividad Aparente (Ohm.m)'] = rho_aparente

        # 3. Graficar los resultados

        fig, ax1 = plt.subplots(figsize=(10, 7))

        # Gráfico principal (Eje Y izquierda: Voltaje)
        color = 'tab:red'
        ax1.set_xlabel('Tiempo (s)', fontsize=12)
        ax1.set_ylabel('Voltaje Inducido (V)', color=color, fontsize=12)
        ax1.loglog(tiempo, np.abs(voltaje), marker='.', linestyle='-', color=color, label='Voltaje Medido')
        ax1.tick_params(axis='y', labelcolor=color)
        ax1.grid(True, which="both", ls="--", alpha=0.5)

        # Escala logarítmica en ambos ejes
        ax1.set_xscale('log')
        ax1.set_yscale('log')
        ax1.invert_xaxis() # Es común invertir el tiempo en TEM, aunque no es estrictamente necesario

        # Eje Y derecha (Resistividad Aparente)
        ax2 = ax1.twinx()
        color = 'tab:blue'
        ax2.set_ylabel('Resistividad Aparente ($\Omega \cdot$m)', color=color, fontsize=12)
        ax2.loglog(tiempo, rho_aparente, marker='o', linestyle='none', color=color, label='Resistividad Aparente')
        ax2.tick_params(axis='y', labelcolor=color)
        ax2.set_yscale('log')

        plt.title('Curva de Decaimiento TEM y Resistividad Aparente', fontsize=14)
        fig.tight_layout()
        plt.show()

        print("\n--- Resultados del Análisis Preliminar ---")
        print(df.head())
        print(f"\nResistividad Aparente Mínima: {rho_aparente.min():.2f} Ohm.m")
        print(f"Resistividad Aparente Máxima: {rho_aparente.max():.2f} Ohm.m")

    except FileNotFoundError:
        print(f"Error: No se encontró el archivo en la ruta especificada: {ruta_archivo}")
    except Exception as e:
        print(f"Ocurrió un error durante el procesamiento: {e}")

# --- Ejemplo de uso ---
# NOTA: Debes reemplazar 'datos_tem.csv' con la ruta real de tu archivo.
# Asegúrate de crear un archivo CSV con las columnas 'Tiempo (s)' y 'Voltaje (V)'.
# analizar_datos_tem('ruta/a/tu/archivo/datos_tem.csv')
