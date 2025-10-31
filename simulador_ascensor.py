#!/usr/bin/env python
import curses
import random
import time
import heapq
from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import List, Optional, Callable, Tuple
from enum import Enum

# Enumeraciones sin cambios
class PersonState(Enum):
    WAITING = "Esperando"
    IN_ELEVATOR = "En ascensor"
    OUT = "Fuera"

class ElevatorDirection(Enum):
    UP = "↑"
    DOWN = "↓"
    IDLE = "Parado"

class ActivityType(Enum):
    WORK = "Trabajo"
    SHOPPING = "Compras"
    EAT = "Comer"
    SLEEP = "Dormir"
    ENTERTAIN = "Entretenerse"
    PARTY = "Fiesta"

# Clase Person sin cambios
@dataclass
class Person:
    """Clase que representa a una persona con un comportamiento aleatorio."""
    id: int
    home_floor: int
    state: PersonState = PersonState.WAITING
    activity: Optional[ActivityType] = None
    activity_time: Optional[float] = None
    destination_floor: Optional[int] = None

    def decide_activity(self, total_floors: int) -> tuple[bool, int, ActivityType]:
        """Decide si la persona sale, a qué piso y para qué actividad."""
        if self.state != PersonState.WAITING:
            return False, -1, None
        rand = random.random()
        if rand < 0.01:  # 1% probabilidad de fiesta
            dest_floor = random.randint(0, total_floors - 1)
            self.activity_time = random.uniform(30, 120)
            return True, dest_floor, ActivityType.PARTY
        elif rand < 0.06:  # 5% probabilidad de trabajo o compras
            self.activity_time = random.uniform(10, 60)
            activity = random.choice([ActivityType.WORK, ActivityType.SHOPPING])
            return True, 0, activity  # Trabajo y compras en piso 0
        return False, -1, None

    def decide_return_activity(self, total_floors: int) -> tuple[bool, int, ActivityType]:
        """Decide la actividad al regresar (comer, dormir, entretenerse)."""
        if self.state != PersonState.OUT or not self.activity_time or time.time() < self.activity_time:
            return False, -1, None
        self.state = PersonState.WAITING
        self.activity_time = random.uniform(10, 60)
        activity = random.choice([ActivityType.EAT, ActivityType.SLEEP, ActivityType.ENTERTAIN])
        return True, self.home_floor, activity

# Clase Elevator corregida
class Elevator:
    """Clase que simula un ascensor con capacidad limitada y movimiento realista."""
    def __init__(self, max_capacity: int = 4, total_floors: int = 5):
        """Inicializa el ascensor con capacidad y número de pisos."""
        self.current_floor = 0
        self.direction = ElevatorDirection.IDLE
        self.passengers: List[Person] = []
        self.max_capacity = max_capacity
        self.total_floors = total_floors
        self.requests: List[Tuple[float, int, Person, int, ActivityType]] = []  # (prioridad, person_id, persona, destino, actividad)

    def add_request(self, person: Person, destination: int, activity: ActivityType, priority: float):
        """Añade una solicitud a la cola de prioridad."""
        heapq.heappush(self.requests, (priority, person.id, person, destination, activity))

    def get_next_request(self) -> Optional[Tuple[Person, int, ActivityType]]:
        """Obtiene la solicitud con mayor prioridad."""
        if not self.requests:
            return None
        _, _, person, destination, activity = heapq.heappop(self.requests)
        return person, destination, activity

    def add_passenger(self, person: Person, destination: int, activity: ActivityType) -> float:
        """Añade un pasajero al ascensor, simulando tiempo de entrada (0.5s)."""
        if len(self.passengers) < self.max_capacity:
            person.state = PersonState.IN_ELEVATOR
            person.destination_floor = destination
            person.activity = activity
            self.passengers.append(person)
            time.sleep(0.5)
            return 0.5
        return 0.0

    def remove_passenger(self, person: Person) -> float:
        """Saca a un pasajero, simulando tiempo de salida (0.5s)."""
        if person in self.passengers:
            self.passengers.remove(person)
            person.state = PersonState.OUT
            person.destination_floor = None
            time.sleep(0.5)
            return 0.5
        return 0.0

    def move_to(self, target_floor: int, update_callback: Callable[['Elevator', List[str]], None]) -> float:
        """Mueve el ascensor al piso objetivo, actualizando la interfaz en cada paso."""
        total_time = 0.0
        new_events = []
        if self.current_floor != target_floor:
            new_events.append(f"Ascensor inicia movimiento desde piso {self.current_floor} a {target_floor}")
            time.sleep(2.0)
            total_time += 2.0
            step = 1 if target_floor > self.current_floor else -1
            while self.current_floor != target_floor:
                self.current_floor += step
                self.direction = ElevatorDirection.UP if step > 0 else ElevatorDirection.DOWN
                new_events.append(f"Ascensor en piso {self.current_floor}")
                update_callback(self, new_events)
                time.sleep(2.0)
                total_time += 2.0
            time.sleep(2.0)
            total_time += 2.0
        self.direction = ElevatorDirection.IDLE
        update_callback(self, new_events)
        return total_time

# Interfaz Display y CursesDisplay (sin cambios significativos, solo ajuste en la visualización de la cola)
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
            waiting = [f"{p.id}({p.activity.value if p.activity else 'Esperando'})" for p in building.floors[floor] if p.state == PersonState.WAITING]
            line = f"Piso {floor}: {waiting}"
            if elevator.current_floor == floor:
                passenger_info = [f"{p.id}({p.activity.value})" for p in elevator.passengers]
                line += f" | Ascensor: {passenger_info} ({elevator.direction.value})"
            self.stdscr.addstr(building.total_floors - floor, 0, line[:width - 1])
        event_start = building.total_floors + 2
        self.stdscr.addstr(event_start, 0, "Eventos:")
        for i, event in enumerate(events[-5:]):
            if event_start + i + 1 < height:
                self.stdscr.addstr(event_start + i + 1, 0, event[:width - 1])
        queue_status = f"Cola: {[(p, d, a.value, prio) for prio, _, p, d, a in elevator.requests]}"
        self.stdscr.addstr(event_start + 6, 0, queue_status[:width - 1])
        self.stdscr.refresh()

# Clase Building corregida
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

    def calculate_priority(self, person: Person, destination: int, activity: ActivityType) -> float:
        """Calcula la prioridad de una solicitud basada en distancia y actividad."""
        distance = abs(self.elevator.current_floor - person.home_floor)
        activity_weights = {
            ActivityType.WORK: 1,
            ActivityType.SHOPPING: 2,
            ActivityType.PARTY: 3,
            ActivityType.EAT: 4,
            ActivityType.SLEEP: 4,
            ActivityType.ENTERTAIN: 4
        }
        activity_priority = activity_weights.get(activity, 4)
        return distance * 10 + activity_priority

    def update_display(self, elevator: Elevator, events: List[str]):
        """Función de devolución para actualizar la interfaz."""
        if self.display:
            self.display.update(self, elevator, self.events + events)

    def run(self):
        """Ejecuta la simulación en un bucle principal."""
        while True:
            self.current_time = time.time()
            new_events = []

            # Generar nuevas solicitudes
            for floor in range(self.total_floors):
                for person in self.floors[floor][:]:
                    if person.state == PersonState.WAITING:
                        should_leave, dest_floor, activity = person.decide_activity(self.total_floors)
                        if should_leave:
                            priority = self.calculate_priority(person, dest_floor, activity)
                            new_events.append(f"Persona {person.id} en piso {floor} solicita ascensor a piso {dest_floor} ({activity.value}, prioridad {priority:.1f})")
                            self.elevator.add_request(person, dest_floor, activity, priority)
                    elif person.state == PersonState.OUT:
                        should_return, dest_floor, activity = person.decide_return_activity(self.total_floors)
                        if should_return:
                            priority = self.calculate_priority(person, dest_floor, activity)
                            new_events.append(f"Persona {person.id} regresa al piso {dest_floor} para {activity.value} (prioridad {priority:.1f})")
                            self.elevator.add_request(person, dest_floor, activity, priority)

            # Procesar solicitudes del ascensor
            request = self.elevator.get_next_request()
            if request:
                person, dest_floor, activity = request
                if person.state == PersonState.WAITING and person in self.floors[person.home_floor]:
                    priority = self.calculate_priority(person, dest_floor, activity)
                    new_events.append(f"Ascensor atiende a persona {person.id} en piso {person.home_floor} ({activity.value}, prioridad {priority:.1f})")
                    self.elevator.move_to(person.home_floor, self.update_display)
                    entry_time = self.elevator.add_passenger(person, dest_floor, activity)
                    if entry_time > 0:
                        new_events.append(f"Persona {person.id} entra al ascensor en piso {person.home_floor} ({activity.value})")
                        self.floors[person.home_floor].remove(person)
                    else:
                        new_events.append(f"Ascensor lleno, persona {person.id} espera")
                        self.elevator.add_request(person, dest_floor, activity, priority)
                    self.elevator.move_to(dest_floor, self.update_display)
                    exit_time = self.elevator.remove_passenger(person)
                    if exit_time > 0:
                        new_events.append(f"Persona {person.id} sale del ascensor en piso {dest_floor} ({activity.value})")
                        self.floors[dest_floor].append(person)
                else:
                    new_events.append(f"Solicitud de persona {person.id} ignorada: estado inválido o no en piso {person.home_floor}")

            # Actualizar eventos
            self.events.extend(new_events)

            # Actualizar pantalla
            self.update_display(self.elevator, [])

            # Pequeña pausa
            time.sleep(0.1)

# Función principal
def main(stdscr):
    building = Building(total_floors=5, display=CursesDisplay(stdscr))
    building.run()

if __name__ == "__main__":
    curses.wrapper(main)
