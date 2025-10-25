import curses
import random
import time
from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import List, Optional
from enum import Enum
import queue

# Enumeraciones sin cambios
class PersonState(Enum):
    WAITING = "Esperando"
    IN_ELEVATOR = "En ascensor"
    OUT = "Fuera"

class ElevatorDirection(Enum):
    UP = "↑"
    DOWN = "↓"
    IDLE = "Parado"

# Clase Person sin cambios
@dataclass
class Person:
    """Clase que representa a una persona con un comportamiento aleatorio."""
    id: int
    home_floor: int
    state: PersonState = PersonState.WAITING
    activity_time: Optional[float] = None

    def decide_activity(self, total_floors: int) -> tuple[bool, int]:
        """Decide si la persona sale y a qué piso (aleatorio, distinto al suyo)."""
        if self.state == PersonState.WAITING and random.random() < 0.05:
            dest_floor = random.randint(0, total_floors - 1)
            while dest_floor == self.home_floor:
                dest_floor = random.randint(0, total_floors - 1)
            self.activity_time = random.uniform(10, 60)
            return True, dest_floor
        return False, -1

    def update(self, current_time: float) -> bool:
        """Actualiza el estado: regresa a casa si termina la actividad."""
        if self.state == PersonState.OUT and self.activity_time:
            if current_time >= self.activity_time:
                self.state = PersonState.WAITING
                self.activity_time = None
                return True
        return False

# Clase Elevator con depuración adicional
class Elevator:
    """Clase que simula un ascensor con capacidad limitada y movimiento realista."""
    def __init__(self, max_capacity: int = 4, total_floors: int = 5):
        """Inicializa el ascensor con capacidad y número de pisos."""
        self.current_floor = 0
        self.direction = ElevatorDirection.IDLE
        self.passengers: List[Person] = []
        self.max_capacity = max_capacity
        self.total_floors = total_floors
        self.requests: queue.Queue = queue.Queue()

    def add_passenger(self, person: Person, destination: int) -> float:
        """Añade un pasajero al ascensor, simulando tiempo de entrada (0.5s)."""
        if len(self.passengers) < self.max_capacity:
            person.state = PersonState.IN_ELEVATOR
            self.passengers.append(person)
            time.sleep(0.5)
            return 0.5
        return 0.0

    def remove_passenger(self, person: Person) -> float:
        """Saca a un pasajero, simulando tiempo de salida (0.5s)."""
        if person in self.passengers:
            self.passengers.remove(person)
            person.state = PersonState.OUT
            time.sleep(0.5)
            return 0.5
        return 0.0

    def move_to(self, target_floor: int) -> float:
        """Mueve el ascensor al piso objetivo, simulando tiempos realistas."""
        total_time = 0.0
        if self.current_floor != target_floor:
            time.sleep(2.0)
            total_time += 2.0
            step = 1 if target_floor > self.current_floor else -1
            while self.current_floor != target_floor:
                time.sleep(2.0)
                total_time += 2.0
                self.current_floor += step
                self.direction = ElevatorDirection.UP if step > 0 else ElevatorDirection.DOWN
            time.sleep(2.0)
            total_time += 2.0
        self.direction = ElevatorDirection.IDLE
        return total_time

# Interfaz Display y CursesDisplay con mejora para mostrar la cola
class Display(ABC):
    """Interfaz para mostrar el estado del edificio y el ascensor."""
    @abstractmethod
    def update(self, building: 'Building', elevator: Elevator, events: List[str]) -> None:
        pass

class CursesDisplay(Display):
    """Muestra el estado del edificio y el ascensor en una interfaz de texto dinámica."""
    def __init__(self, stdscr):
        """Inicializa la ventana de curses."""
        self.stdscr = stdscr
        curses.curs_set(0)
        self.stdscr.timeout(100)

    def update(self, building: 'Building', elevator: Elevator, events: List[str]) -> None:
        """Actualiza la pantalla con el estado del edificio, ascensor y eventos."""
        self.stdscr.clear()
        height, width = self.stdscr.getmaxyx()
        for floor in range(building.total_floors - 1, -1, -1):
            waiting = [p.id for p in building.floors[floor] if p.state == PersonState.WAITING]
            line = f"Piso {floor}: {waiting}"
            if elevator.current_floor == floor:
                passenger_ids = [p.id for p in elevator.passengers]
                line += f" | Ascensor: {passenger_ids} ({elevator.direction.value})"
            self.stdscr.addstr(building.total_floors - floor, 0, line[:width - 1])
        event_start = building.total_floors + 2
        self.stdscr.addstr(event_start, 0, "Eventos:")
        for i, event in enumerate(events[-5:]):
            if event_start + i + 1 < height:
                self.stdscr.addstr(event_start + i + 1, 0, event[:width - 1])
        # Mostrar estado de la cola
        queue_status = f"Cola: {list(elevator.requests.queue)}"
        self.stdscr.addstr(event_start + 6, 0, queue_status[:width - 1])
        self.stdscr.refresh()

# Clase Building corregida para procesar solicitudes correctamente
class Building:
    """Clase que representa el edificio y gestiona la simulación."""
    def __init__(self, total_floors: int = 5, display: Display = None):
        """Inicializa el edificio con pisos, personas, ascensor y visualización."""
        self.total_floors = total_floors
        self.floors: List[List[Person]] = [[] for _ in range(total_floors)]
        self.elevator = Elevator(max_capacity=4, total_floors=total_floors)
        self.display = display
        self.events: List[str] = []
        self.current_time = time.time()
        self.person_id = 0
        for floor in range(total_floors):
            num_people = random.randint(1, 5)
            for _ in range(num_people):
                self.floors[floor].append(Person(id=self.person_id, home_floor=floor))
                self.person_id += 1

    def run(self):
        """Ejecuta la simulación en un bucle principal."""
        while True:
            self.current_time = time.time()
            new_events = []

            # Generar nuevas solicitudes
            for floor in range(self.total_floors):
                for person in self.floors[floor][:]:
                    if person.state == PersonState.WAITING:
                        should_leave, dest_floor = person.decide_activity(self.total_floors)
                        if should_leave:
                            new_events.append(f"Persona {person.id} en piso {floor} solicita ascensor a piso {dest_floor}")
                            self.elevator.requests.put((person, dest_floor))
                    elif person.state == PersonState.OUT:
                        if person.update(self.current_time):
                            new_events.append(f"Persona {person.id} regresa al piso {person.home_floor}")
                            self.elevator.requests.put((person, person.home_floor))

            # Procesar solicitudes del ascensor
            if not self.elevator.requests.empty():
                person, dest_floor = self.elevator.requests.get()
                # Verificar que la solicitud es válida
                if person.state == PersonState.WAITING and person in self.floors[person.home_floor]:
                    new_events.append(f"Ascensor se mueve al piso {person.home_floor}")
                    self.elevator.move_to(person.home_floor)
                    # Persona entra al ascensor
                    entry_time = self.elevator.add_passenger(person, dest_floor)
                    if entry_time > 0:
                        new_events.append(f"Persona {person.id} entra al ascensor en piso {person.home_floor}")
                        self.floors[person.home_floor].remove(person)
                    else:
                        new_events.append(f"Ascensor lleno, persona {person.id} espera")
                        self.elevator.requests.put((person, dest_floor))  # Reinsertar solicitud
                    # Mover al destino
                    new_events.append(f"Ascensor se mueve al piso {dest_floor}")
                    self.elevator.move_to(dest_floor)
                    # Persona sale del ascensor
                    exit_time = self.elevator.remove_passenger(person)
                    if exit_time > 0:
                        new_events.append(f"Persona {person.id} sale del ascensor en piso {dest_floor}")
                        self.floors[dest_floor].append(person)
                else:
                    new_events.append(f"Solicitud de persona {person.id} ignorada: estado inválido o no en piso {person.home_floor}")

            # Actualizar eventos
            self.events.extend(new_events)

            # Actualizar pantalla
            if self.display:
                self.display.update(self, self.elevator, self.events)

            # Pequeña pausa para evitar uso excesivo de CPU
            time.sleep(0.1)

# Función principal sin cambios
def main(stdscr):
    building = Building(total_floors=5, display=CursesDisplay(stdscr))
    building.run()

if __name__ == "__main__":
    curses.wrapper(main)
