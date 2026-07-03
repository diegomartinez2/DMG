import datetime
import tkinter as tk
from tkinter import messagebox


# ==========================================
# 1. CAPA DE LÓGICA (K.I.S.S. & Datos)
# ==========================================
class RegistroManager:
    """Gestiona el estado interno de ocupación y los registros de tiempo."""

    def __init__(self, capacidad_maxima: int):
        self.capacidad_maxima = capacidad_maxima
        # Diccionario para rastrear quién ocupa qué número de slot {slot_id: "Nombre"}
        self.ocupantes = {}

    def registrar_entrada(self, slot_id: int, nombre: str) -> str:
        """Guarda la hora de entrada y ocupa el slot."""
        nombre_limpio = nombre.strip()
        if not nombre_limpio:
            raise ValueError("El nombre no puede estar vacío.")

        ahora = datetime.datetime.now().strftime("%H:%M:%S")
        self.ocupantes[slot_id] = nombre_limpio

        # Aquí podrías añadir código para guardar en un archivo TXT o Base de datos
        print(f"[ENTRADA] {nombre_limpio} - Slot {slot_id} - Hora: {ahora}")
        return ahora

    def registrar_salida(self, slot_id: int) -> str:
        """Guarda la hora de salida y libera el slot."""
        nombre = self.ocupantes.get(slot_id, "Desconocido")
        ahora = datetime.datetime.now().strftime("%H:%M:%S")

        print(f"[SALIDA] {nombre} - Slot {slot_id} - Hora: {ahora}")

        # Liberamos el slot del registro activo
        if slot_id in self.ocupantes:
            del self.ocupantes[slot_id]

        return ahora

    @property
    def personas_actuales(self) -> int:
        return len(self.ocupantes)


# ==========================================
# 2. CAPA DE INTERFAZ MODULAR (Componentes)
# ==========================================
class FilaPersonaUI(tk.Frame):
    """Componente visual reutilizable para cada fila de registro."""

    def __init__(self, parent, slot_id: int, manager: RegistroManager, *args, **kwargs):
        super().__init__(parent, *args, **kwargs)
        self.slot_id = slot_id
        self.manager = manager

        # Configuración de diseño de la fila
        self.columnconfigure(1, weight=1)

        # 1. Indicador de número de persona
        self.lbl_numero = tk.Label(
            self, text=f"Persona {slot_id}:", width=10, anchor="w"
        )
        self.lbl_numero.grid(row=0, column=0, padx=5, pady=5)

        # 2. Campo de texto para el nombre
        self.ent_nombre = tk.Entry(self, font=("Arial", 10))
        self.ent_nombre.grid(row=0, column=1, padx=5, pady=5, sticky="ew")

        # 3. Botón de Entrada
        self.btn_entrada = tk.Button(
            self,
            text="Entrada",
            bg="#4CAF50",
            fg="white",
            width=8,
            command=self._ejecutar_entrada,
        )
        self.btn_entrada.grid(row=0, column=2, padx=5, pady=5)

        # 4. Botón de Salida (Empieza desactivado)
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

        # 5. Etiqueta de estado de tiempos
        self.lbl_tiempos = tk.Label(
            self, text="Esperando entrada...", fg="gray", width=25, anchor="w"
        )
        self.lbl_tiempos.grid(row=0, column=4, padx=10, pady=5)

    def _ejecutar_entrada(self):
        nombre = self.ent_nombre.get()
        try:
            hora = self.manager.registrar_entrada(self.slot_id, nombre)

            # Modificaciones de interfaz tras el éxito
            self.lbl_tiempos.config(text=f"Entrada: {hora}", fg="green")
            self.ent_nombre.config(state=tk.DISABLED)  # Bloquea el nombre
            self.btn_entrada.config(state=tk.DISABLED)  # Bloquea botón entrada
            self.btn_salida.config(state=tk.NORMAL)  # Activa botón salida
        except ValueError as e:
            messagebox.showwarning("Error de Entrada", str(e))

    def _ejecutar_salida(self):
        hora_entrada = self.lbl_tiempos.cget("text").split(": ")[1]
        hora_salida = self.manager.registrar_salida(self.slot_id)

        # Modificaciones de interfaz tras la salida
        self.lbl_tiempos.config(
            text=f"In: {hora_entrada} | Out: {hora_salida}", fg="blue"
        )
        self.btn_salida.config(state=tk.DISABLED)

        # Pequeño reset para que el slot pueda volver a usarse si se vacía la fila
        self.ent_nombre.config(state=tk.NORMAL)
        self.ent_nombre.delete(0, tk.END)
        self.btn_entrada.config(state=tk.NORMAL)


# ==========================================
# 3. APLICACIÓN PRINCIPAL (Orquestador)
# ==========================================
class AplicacionSala:

    def __init__(self, root, capacidad_sala: int = 6):
        self.root = root
        self.root.title("Control de Aforo y Registro de Sala")
        self.root.geometry("650x350")
        self.root.resizable(False, False)

        # Inicializamos la lógica con la capacidad deseada
        self.manager = RegistroManager(capacidad_maxima=capacidad_sala)

        # Título superior
        lbl_titulo = tk.Label(
            root,
            text=f"Registro de Tiempos (Capacidad Máxima: {capacidad_sala} personas)",
            font=("Arial", 12, "bold"),
            pady=10,
        )
        lbl_titulo.pack()

        # Contenedor para las filas modulares
        self.frame_filas = tk.Frame(root)
        self.frame_filas.pack(padx=15, pady=5, fill="x", expand=True)

        # Bucle dinámico: Genera tantas filas como indique 'capacidad_sala'
        self.filas = []
        for i in range(1, capacidad_sala + 1):
            fila = FilaPersonaUI(self.frame_filas, slot_id=i, manager=self.manager)
            fila.pack(fill="x", pady=2)
            self.filas.append(fila)


# ==========================================
# EJECUCIÓN DEL PROGRAMA
# ==========================================
if __name__ == "__main__":
    # CONFIGURACIÓN FUTURA: Si el día de mañana cambias este 6 por un 10 o un 4,
    # la interfaz se adaptará y crecerá automáticamente sin tocar nada más.
    CAPACIDAD_CONFIGURADA = 6

    root = tk.Tk()
    app = AplicacionSala(root, capacidad_sala=CAPACIDAD_CONFIGURADA)
    root.mainloop()
