import random
import time
from typing import Dict, List, Optional, Tuple

# Importaciones para la GUI y la gráfica
import tkinter as tk
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.cm as cm
import matplotlib.pyplot as plt # Necesario para plt.Line2D o Line2D

# --- Parámetros de la Simulación ---
TASA_BASE_REPLICACION = 0.04  # Probabilidad base de que un individuo se replique
TASA_BASE_MORTALIDAD = 0.02   # Probabilidad base de que un individuo muera
TASA_BASE_MUTACION = 0.04     # Probabilidad base de que una réplica sea una mutación
TASA_CREACION_ESPONTANEA = 0.01 # SOLO para el Elemento 0

# Configuración del entorno
NUMERO_CICLOS = 400            # Aumentado para mejor visualización
POBLACION_INICIAL_ELEMENTO_0 = 100
MAX_DESVIACION_MUTACION = 0.10
MAX_POBLACION_TOTAL = 1000000
INTERVALO_MS = 100             # Retraso en milisegundos entre cada ciclo de actualización de la gráfica

# --- Clase Elemento: Define el Tipo de Organismo ---

class Elemento:
    """
    Define un tipo de elemento (o especie) con sus tasas de vida y evolución.
    """
    next_id: int = 0

    def __init__(self,
                 tasa_replicacion: float,
                 tasa_mortalidad: float,
                 tasa_mutacion: float,
                 generacion: int,
                 parent_id: Optional[int] = None):
        """Inicializa un nuevo tipo de elemento."""
        self.id = Elemento.next_id
        Elemento.next_id += 1

        self.tasa_replicacion = tasa_replicacion
        self.tasa_mortalidad = tasa_mortalidad
        self.tasa_mutacion = tasa_mutacion
        self.generacion = generacion
        self.parent_id = parent_id

    def __repr__(self):
        """Representación legible para el seguimiento."""
        parent_info = f"Padre: {self.parent_id}" if self.parent_id is not None else "Origen"
        return (f"E{self.id} (G{self.generacion}, {parent_info}): "
                f"R={self.tasa_replicacion:.3f}, M={self.tasa_mortalidad:.3f}")

# --- Clase Simulación: Controla el Ecosistema ---

class Simulacion:
    """
    Gestiona los ciclos de vida y la evolución de los elementos.
    Refactorizada para ejecutar un ciclo a la vez.
    """

    def __init__(self, poblacion_inicial: int):
        self.tipos_elementos: Dict[int, Elemento] = {}
        self.poblacion_actual: Dict[int, int] = {}
        self.ciclo_actual: int = 0
        self.total_poblacion_historico: List[int] = [] # Para la gráfica principal
        self.poblaciones_por_tipo: Dict[int, List[int]] = {} # Para la gráfica de especies

        # Crear el Elemento 0
        self.elemento_base = Elemento(
            tasa_replicacion=TASA_BASE_REPLICACION,
            tasa_mortalidad=TASA_BASE_MORTALIDAD,
            tasa_mutacion=TASA_BASE_MUTACION,
            generacion=0
        )
        self.tipos_elementos[self.elemento_base.id] = self.elemento_base
        self.poblacion_actual[self.elemento_base.id] = poblacion_inicial

        # El historial del Elemento 0 incluye su población inicial (ciclo 0)
        self.poblaciones_por_tipo[self.elemento_base.id] = [poblacion_inicial]
        self.total_poblacion_historico.append(poblacion_inicial)

    def _mutar_tasa(self, tasa_original: float) -> float:
        """Calcula una nueva tasa mutada, ligeramente aleatoria."""
        desviacion = random.uniform(-MAX_DESVIACION_MUTACION, MAX_DESVIACION_MUTACION)
        nueva_tasa = tasa_original * (1.0 + desviacion)
        return max(0.001, min(1.0, nueva_tasa))

    def _crear_mutacion(self, padre_id: int, elemento_padre: Elemento) -> Elemento:
        """Genera un nuevo tipo de elemento mutado."""
        nueva_tasa_rep = self._mutar_tasa(elemento_padre.tasa_replicacion)
        nueva_tasa_mort = self._mutar_tasa(elemento_padre.tasa_mortalidad)
        nueva_tasa_mut = self._mutar_tasa(elemento_padre.tasa_mutacion)

        nuevo_elemento = Elemento(
            tasa_replicacion=nueva_tasa_rep,
            tasa_mortalidad=nueva_tasa_mort,
            tasa_mutacion=nueva_tasa_mut,
            generacion=elemento_padre.generacion + 1,
            parent_id=padre_id
        )

        self.tipos_elementos[nuevo_elemento.id] = nuevo_elemento
        self.poblacion_actual[nuevo_elemento.id] = 1

        # INICIALIZACIÓN CORRECTA:
        # La lista debe tener ceros para los ciclos 0 hasta (ciclo_actual - 1),
        # y luego 1 para el ciclo_actual.
        # self.ciclo_actual ya fue incrementado en ejecutar_ciclo antes de esta llamada.
        # Por ejemplo: Si muta en ciclo 5, la lista debe ser [0, 0, 0, 0, 0, 1] (longitud 6)
        self.poblaciones_por_tipo[nuevo_elemento.id] = [0] * self.ciclo_actual + [1]

        return nuevo_elemento

    def _procesar_poblacion(self, id_tipo: int, poblacion: int, factor_densidad: float) -> Dict[int, int]:
        """Calcula los cambios (muertes, réplicas, mutaciones) para un tipo de elemento."""
        elemento = self.tipos_elementos[id_tipo]
        cambios = {id_tipo: 0}

        tasa_replicacion_efectiva = elemento.tasa_replicacion * factor_densidad

        for _ in range(poblacion):
            if random.random() < elemento.tasa_mortalidad:
                cambios[id_tipo] -= 1
                continue

            if random.random() < tasa_replicacion_efectiva:
                if random.random() < elemento.tasa_mutacion:
                    self._crear_mutacion(id_tipo, elemento)
                else:
                    cambios[id_tipo] += 1

        return cambios

    def _creacion_espontanea(self, factor_densidad: float) -> Dict[int, int]:
        """Añade individuos del Elemento 0, afectado por la limitación de capacidad."""
        cambio_espontaneo = 0
        tasa_efectiva = TASA_CREACION_ESPONTANEA * factor_densidad

        if random.random() < tasa_efectiva:
            cambio_espontaneo = 1

        return {0: cambio_espontaneo}

    def _aplicar_cambios(self, cambios: Dict[int, int]):
        """Aplica todos los cambios calculados a la población actual."""
        for id_tipo, delta in cambios.items():
            if id_tipo not in self.poblacion_actual:
                self.poblacion_actual[id_tipo] = 0

            self.poblacion_actual[id_tipo] += delta

            if self.poblacion_actual[id_tipo] <= 0:
                # El elemento se extingue
                del self.poblacion_actual[id_tipo]


    def ejecutar_ciclo(self) -> Tuple[bool, int, str]:
        """
        Ejecuta un único ciclo de la simulación.
        Retorna (continuar_simulacion, población_total_anterior, mensaje_de_log).
        """
        # 1. Preparación
        self.ciclo_actual += 1
        cambios_totales = {}
        poblacion_anterior = sum(self.poblacion_actual.values())
        mensaje = ""

        if self.ciclo_actual > NUMERO_CICLOS or poblacion_anterior == 0:
            return (False, poblacion_anterior, f"Simulación finalizada en ciclo {self.ciclo_actual - 1} o Extinción total.")

        # 2. Factor de Densidad
        factor_densidad = max(0.0, 1.0 - (poblacion_anterior / MAX_POBLACION_TOTAL))

        if factor_densidad == 0.0:
            mensaje += f"Ciclo {self.ciclo_actual}: Capacidad Límite (K) alcanzada.\n"

        # 3. Procesar Cambios
        cambios_esp = self._creacion_espontanea(factor_densidad)
        cambios_totales.update(cambios_esp)

        for id_tipo, poblacion in list(self.poblacion_actual.items()):
            cambios_vida = self._procesar_poblacion(id_tipo, poblacion, factor_densidad)

            for k, v in cambios_vida.items():
                cambios_totales[k] = cambios_totales.get(k, 0) + v

        # 4. Aplicar Cambios
        self._aplicar_cambios(cambios_totales)

        # 5. Registro de Resultados
        total_poblacion = sum(self.poblacion_actual.values())
        diferencia = total_poblacion - poblacion_anterior

        nuevos_tipos = [e for e in self.tipos_elementos.values() if e.id in self.poblacion_actual and e.generacion == self.ciclo_actual]

        if nuevos_tipos:
            mensaje += f"¡Mutación(es) en G{self.ciclo_actual}!\n"
            for nt in nuevos_tipos:
                 mensaje += f"  + E{nt.id} (R={nt.tasa_replicacion:.3f}, M={nt.tasa_mortalidad:.3f})\n"

        mensaje += f"Ciclo {self.ciclo_actual:03d} | Pop Total: {total_poblacion:,} | Cambio: {diferencia:+}"

        # 6. Actualizar historial de poblaciones (para la gráfica)
        self.total_poblacion_historico.append(total_poblacion)

        tipos_existentes = set(self.tipos_elementos.keys())

        for tipo_id in tipos_existentes:
            if tipo_id not in self.poblaciones_por_tipo:
                # Esto no debería pasar con la lógica _crear_mutacion, pero es un buen control
                continue

            # Obtener la población actual o 0 si se extinguió
            pop_actual_tipo = self.poblacion_actual.get(tipo_id, 0)

            # Añadir la población del ciclo actual (Efecto Inmediato)
            # Ya sea que el tipo haya existido siempre, haya mutado en este ciclo, o se haya extinguido
            self.poblaciones_por_tipo[tipo_id].append(pop_actual_tipo)

        return (True, poblacion_anterior, mensaje)

# --- Clase de la Aplicación con GUI y Gráfica ---

class AppGrafica:
    """
    Clase principal que maneja la interfaz gráfica Tkinter y la visualización
    de la simulación de la vida estocástica.
    """
    def __init__(self, master):
        self.master = master
        master.title("Simulación de Vida con Gráfica Dinámica (Logística)")
        self.simulacion = Simulacion(POBLACION_INICIAL_ELEMENTO_0)
        self.corriendo = False
        self.ciclo_actual = 0

        # Configuración de la gráfica (Actualizado para compatibilidad con Matplotlib 3.8+)
        self.fig = Figure(figsize=(8, 6), dpi=100)
        self.ax = self.fig.add_subplot(111)
        self.canvas = FigureCanvasTkAgg(self.fig, master=master)
        self.canvas_widget = self.canvas.get_tk_widget()
        self.canvas_widget.pack(side=tk.TOP, fill=tk.BOTH, expand=1)

        # Usar matplotlib.colormaps o pyplot.get_cmap como se sugiere en el warning
        self.colores = plt.colormaps.get_cmap('hsv', 20) # Mapa de colores para las especies
        self.lineas: Dict[int, plt.Line2D] = {} # Almacena las líneas de la gráfica

        # Inicialización de la UI
        self._inicializar_ui()
        self._inicializar_grafica()

    def _inicializar_ui(self):
        """Inicializa los elementos de la interfaz de usuario (Tkinter)."""
        frame = tk.Frame(self.master)
        frame.pack(side=tk.BOTTOM, fill=tk.X)

        # Etiqueta de estado
        self.estado_var = tk.StringVar(self.master, value="Lista para Iniciar")
        self.estado_label = tk.Label(frame, textvariable=self.estado_var, padx=10, pady=5)
        self.estado_label.pack(side=tk.LEFT)

        # Botón de inicio/pausa
        self.start_button = tk.Button(frame, text="Iniciar Simulación", command=self.toggle_simulacion)
        self.start_button.pack(side=tk.LEFT, padx=5)

        # Área de log/información
        self.log_text = tk.Text(self.master, height=10, state=tk.DISABLED)
        self.log_text.pack(side=tk.BOTTOM, fill=tk.X, padx=5, pady=5)
        self.log_mensaje(f"Límite de Capacidad (K): {MAX_POBLACION_TOTAL:,} individuos.")
        self.log_mensaje(f"Inicia E0: {POBLACION_INICIAL_ELEMENTO_0} individuos.")

    def _inicializar_grafica(self):
        """Configura los ejes y títulos de la gráfica."""
        self.ax.set_title("Crecimiento Poblacional Estocástico (N vs K)")
        self.ax.set_xlabel("Ciclo de Tiempo")
        self.ax.set_ylabel("Población")
        self.ax.axhline(MAX_POBLACION_TOTAL, color='red', linestyle='--', label='Capacidad Límite (K)')
        self.ax.legend(loc='upper left')
        self.ax.set_xlim(0, NUMERO_CICLOS)
        self.ax.set_ylim(0, MAX_POBLACION_TOTAL * 1.1)
        self.fig.tight_layout()

    def log_mensaje(self, mensaje: str):
        """Añade mensajes al área de texto del log."""
        self.log_text.config(state=tk.NORMAL)
        self.log_text.insert(tk.END, mensaje + "\n")
        self.log_text.see(tk.END)
        self.log_text.config(state=tk.DISABLED)

    def toggle_simulacion(self):
        """Alterna el estado de la simulación (corriendo/pausada)."""
        self.corriendo = not self.corriendo
        if self.corriendo:
            self.start_button.config(text="Pausar Simulación")
            # Iniciar la primera ejecución después de un breve retraso
            self.master.after(INTERVALO_MS, self.ejecutar_y_actualizar)
        else:
            self.start_button.config(text="Reanudar Simulación")

    def ejecutar_y_actualizar(self):
        """Ejecuta un ciclo de la simulación y actualiza la gráfica."""
        if not self.corriendo:
            return

        continuar, pop_anterior, mensaje_log = self.simulacion.ejecutar_ciclo()
        self.ciclo_actual = self.simulacion.ciclo_actual
        self.log_mensaje(mensaje_log)
        self.estado_var.set(f"Ciclo: {self.ciclo_actual}/{NUMERO_CICLOS} | Total: {sum(self.simulacion.poblacion_actual.values()):,}")

        # La lista de ciclos debe tener longitud (ciclo_actual + 1)
        # Si ciclo_actual es 5, queremos ciclos [0, 1, 2, 3, 4, 5] (longitud 6)
        ciclos = list(range(self.ciclo_actual + 1))

        # Crear o actualizar líneas
        tipos_a_mostrar = list(self.simulacion.tipos_elementos.keys())

        # Iterar sobre todos los tipos de elementos que han existido
        for tipo_id in tipos_a_mostrar:

            pops = self.simulacion.poblaciones_por_tipo.get(tipo_id, [])

            # --- CORRECCIÓN CRÍTICA ---
            # Asegurar que pops_data tiene la misma longitud que ciclos
            longitud_necesaria = len(ciclos)
            pops_data = pops + [0] * (longitud_necesaria - len(pops))
            # --------------------------

            if tipo_id not in self.lineas:
                # Nuevo tipo de elemento (Mutación)
                color = self.colores(tipo_id % 20)
                elemento = self.simulacion.tipos_elementos[tipo_id]
                etiqueta = f"E{tipo_id} (G{elemento.generacion})"

                # Usar plt.Line2D directamente como tipo de objeto
                linea, = self.ax.plot(ciclos, pops_data, label=etiqueta, color=color, linewidth=1.5, alpha=0.8)
                self.lineas[tipo_id] = linea
                # Forzar la actualización de la leyenda para que aparezca el nuevo elemento
                self.ax.legend(loc='upper left')
            else:
                # Tipo ya existente, solo actualizar datos
                self.lineas[tipo_id].set_data(ciclos, pops_data)

        # 2. Re-dibujar
        self.ax.relim()
        self.ax.autoscale_view()
        self.ax.set_xlim(0, NUMERO_CICLOS) # Mantenemos el límite X fijo
        self.ax.set_ylim(0, MAX_POBLACION_TOTAL * 1.1)
        self.canvas.draw_idle()

        # 3. Planificar la próxima ejecución
        if continuar:
            self.master.after(INTERVALO_MS, self.ejecutar_y_actualizar)
        else:
            self.start_button.config(text="Simulación Finalizada", state=tk.DISABLED)
            self.log_mensaje("\n--- Resumen Final ---")

# --- Punto de Entrada de la Aplicación ---
if __name__ == "__main__":
    # Inicializar la aplicación Tkinter
    root = tk.Tk()

    # Crear y correr la aplicación
    app = AppGrafica(root)
    root.mainloop()
