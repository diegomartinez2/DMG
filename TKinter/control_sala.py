import csv
import datetime
import os
import tkinter as tk
from tkinter import messagebox


# ==========================================
# 1. CAPA DE LÓGICA (Gestión de Datos y CSV)
# ==========================================
class RegistroManager:
    """Gestiona el estado interno de ocupación y escribe en el archivo CSV."""

    def __init__(self, capacidad_maxima: int, nombre_archivo: str = "registro_sala.csv"):
        self.capacidad_maxima = capacidad_maxima
        self.nombre_archivo = nombre_archivo
        self.ocupantes = {}  # Rastrear ocupantes activos {slot_id: "Nombre"}
        self._inicializar_csv()

    def _inicializar_csv(self):
        """Crea el archivo CSV con sus cabeceras si no existe en el sistema."""
        if not os.path.exists(self.nombre_archivo):
            with open(
                self.nombre_archivo, mode="w", newline="", encoding="utf-8"
            ) as f:
                escritor = csv.writer(f)
                escritor.writerow(["Fecha", "Hora", "Accion", "Slot", "Nombre"])

    def _guardar_en_csv(self, accion: str, slot_id: int, nombre: str):
        """Escribe una nueva fila en el registro histórico CSV."""
        ahora = datetime.datetime.now()
        fecha = ahora.strftime("%Y-%m-%d")
        hora = ahora.strftime("%H:%M:%S")

        with open(self.nombre_archivo, mode="a", newline="", encoding="utf-8") as f:
            escritor = csv.writer(f)
            escritor.writerow([fecha, hora, accion, slot_id, nombre])

        return hora

    def registrar_entrada(self, slot_id: int, nombre: str) -> str:
        """Registra la entrada en memoria y en el CSV."""
        nombre_limpio = nombre.strip()
        if not nombre_limpio:
            raise ValueError("El nombre no puede estar vacío.")

        # Guardamos en el diccionario de la sala
        self.ocupantes[slot_id] = nombre_limpio

        # Guardamos en el archivo físico
        hora_entrada = self._guardar_en_csv("ENTRADA", slot_id, nombre_limpio)
        return hora_entrada

    def registrar_salida(self, slot_id: int) -> str:
        """Registra la salida en el CSV y libera el slot de la memoria."""
        nombre = self.ocupantes.get(slot_id, "Desconocido")

        # Guardamos en el archivo físico
        hora_salida = self._guardar_en_csv("SALIDA", slot_id, nombre)

        # Liberamos el slot para que pueda volver a usarse
        if slot_id in self.ocupantes:
            del self.ocupantes[slot_id]

        return hora_salida


# ==========================================
# 2. CAPA DE INTERFAZ MODULAR (Componentes)
# ==========================================
class FilaPersonaUI(tk.Frame):
    """Componente visual reutilizable para cada fila de registro."""

    def __init__(self, parent, slot_id: int, manager: RegistroManager, *args, **kwargs):
        super().__init__(parent, *args, **kwargs)
        self.slot_id = slot_id
        self.manager = manager

        self.columnconfigure(1, weight=1)

        # UI Elements
        self.lbl_numero = tk.Label(
            self, text=f"Persona {slot_id}:", width=10, anchor="w"
        )
        self.lbl_numero.grid(row=0, column=0, padx=5, pady=5)

        self.ent_nombre = tk.Entry(self, font=("Arial", 10))
        self.ent_nombre.grid(row=0, column=1, padx=5, pady=5, sticky="ew")

        self.btn_entrada = tk.Button(
            self,
            text="Entrada",
            bg="#4CAF50",
            fg="white",
            width=8,
            command=self._ejecutar_entrada,
        )
        self.btn_entrada.grid(row=0, column=2, padx=5, pady=5)

        self.btn_salida = tk.Button(
            self,
            text="Salida",
            bg="#f44336",
            fg="white",
            width=8,
            state=tk.DISABLED,
            command=self._ejecutar_salida,
        )
        self.btn_salida.grid(row=0, column=3, padx=5, pady=5)

        self.lbl_tiempos = tk.Label(
            self, text="Esperando entrada...", fg="gray", width=25, anchor="w"
        )
        self.lbl_tiempos.grid(row=0, column=4, padx=10, pady=5)

    def _ejecutar_entrada(self):
        nombre = self.ent_nombre.get()
        try:
            hora = self.manager.registrar_entrada(self.slot_id, nombre)
            self.lbl_tiempos.config(text=f"Entrada: {hora}", fg="green")
            self.ent_nombre.config(state=tk.DISABLED)
            self.btn_entrada.config(state=tk.DISABLED)
            self.btn_salida.config(state=tk.NORMAL)
        except ValueError as e:
            messagebox.showwarning("Error de Entrada", str(e))

    def _ejecutar_salida(self):
        hora_entrada = self.lbl_tiempos.cget("text").split(": ")[1]
        hora_salida = self.manager.registrar_salida(self.slot_id)

        self.lbl_tiempos.config(
            text=f"In: {hora_entrada} | Out: {hora_salida}", fg="blue"
        )
        self.btn_salida.config(state=tk.DISABLED)

        # Reset para permitir que otro usuario ocupe este slot
        self.ent_nombre.config(state=tk.NORMAL)
        self.ent_nombre.delete(0, tk.END)
        self.btn_entrada.config(state=tk.NORMAL)


# ==========================================
# 3. APLICACIÓN PRINCIPAL (Orquestador)
# ==========================================
class AplicacionSala:

    def __init__(self, root, capacidad_sala: int = 6):
        self.root = root
        self.root.title("Control de Aforo e Historial CSV")
        self.root.geometry("650x350")
        self.root.resizable(False, False)

        self.manager = RegistroManager(capacidad_maxima=capacidad_sala)

        lbl_titulo = tk.Label(
            root,
            text=f"Registro de Tiempos (Capacidad Máxima: {capacidad_sala} personas)",
            font=("Arial", 12, "bold"),
            pady=10,
        )
        lbl_titulo.pack()

        self.frame_filas = tk.Frame(root)
        self.frame_filas.pack(padx=15, pady=5, fill="x", expand=True)

        self.filas = []
        for i in range(1, capacidad_sala + 1):
            fila = FilaPersonaUI(self.frame_filas, slot_id=i, manager=self.manager)
            fila.pack(fill="x", pady=2)
            self.filas.append(fila)


# ==========================================
# EJECUCIÓN DEL PROGRAMA
# ==========================================
if __name__ == "__main__":
    # Ajusta este número según las necesidades de la sala
    CAPACIDAD_CONFIGURADA = 6

    root = tk.Tk()
    app = AplicacionSala(root, capacidad_sala=CAPACIDAD_CONFIGURADA)
    root.mainloop()
