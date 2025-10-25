import curses
import random
import time
import heapq
from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import List, Optional, Callable, Tuple
from enum import Enum

# Enumeraciones
class PersonState(Enum):
    WAITING = "Esperando"
    IN_ELEVATOR = "En ascensor"
    OUT = "Fuera"

class ElevatorDirection(Enum):
    UP = "↑"
    DOWN = "↓"
    IDLE = "Parado"

class ActivityType(Enum):
    GOING_TO_WORK = "Yendo a trabajar"
    WORKING = "Trabajando"
    SHOPPING = "Compras"
    EAT = "Comer"
    SLEEP = "Dormir"
    ENTERTAIN = "Entretenerse"
    PARTY = "Fiesta"

# Clase Person
@dataclass
class Person:
    """Clase que representa a una persona con un comportamiento aleatorio."""
    id: int
    home_floor: int
    state: PersonState = PersonState.WAITING
    activity: Optional[ActivityType] = None
    activity_time: Optional[float] = None
    destination_floor: Optional[int] = None
    request_time: Optional[float] = None  # Tiempo en que se generó la solicitud

    def decide_activity(self, total_floors: int) -> tuple[bool, int, ActivityType]:
        """Decide si la persona sale, a qué piso y para qué actividad."""
        if self.state != PersonState.WAITING:
            return False, -1, None
        rand = random.random()
        self.request_time = time.time()
        if rand < 0.01:  # 1% probabilidad de fiesta
            dest_floor = random.randint(0, total_floors - 1)
            self.activity_time = random.uniform(30, 120)
            return True, dest_floor, ActivityType.PARTY
        elif rand < 0.06:  # 5% probabilidad de trabajo o compras
            self.activity_time = random.uniform(10, 60)
            activity = random.choice([ActivityType.GOING_TO_WORK, ActivityType.SHOPPING])
            return True, 0, activity  # Trabajo y compras en piso 0
        elif rand < 0.11:  # 5% probabilidad de entretenimiento
            dest_floor = random.choice([0] + [self.home_floor])  # Piso 0 o piso de origen
            self.activity_time = random.uniform(10, 60)
            return True, dest_floor, ActivityType.ENTERTAIN
        return False, -1, None

    def decide_return_activity(self, total_floors: int) -> tuple[bool, int, ActivityType]:
        """Decide la actividad al regresar (comer, dormir, no en piso 0)."""
        if self.state != PersonState.OUT or not self.activity_time or time.time() < self.activity_time:
            return False, -1, None
        self.state = PersonState.WAITING
        self.request_time = time.time()
        self.activity_time = random.uniform(10, 60)
        activity = random.choice([ActivityType.EAT, ActivityType.SLEEP])  # No ENTERTAIN en retorno
        return True, self.home_floor, activity  # Siempre regresar al piso de origen

    def get_waiting_display(self) -> str:
        """Muestra la actividad deseada con puntos suspensivos según el tiempo de espera."""
        if not self.activity or not self.request_time:
            return "Esperando"
        wait_time = time.time() - self.request_time
        dots = "." * (1 + int(wait_time // 5))  # Un punto más por cada 5s
        return f"{self.activity.value}{dots}"

# Clase Elevator
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

    def get_compatible_requests(self, current_direction: Optional[ElevatorDirection], prioritize_work: bool, floors: List[List[Person]]) -> List[Tuple[Person, int, ActivityType]]:
        """Obtiene solicitudes compatibles con la dirección actual y prioridad de trabajo."""
        compatible_requests = []
        work_requests = []
        for priority, _, person, destination, activity in self.requests[:]:
            if person.state != PersonState.WAITING or person not in floors[person.home_floor]:
                continue
            # Verificar si el destino está en la dirección del movimiento
            if current_direction == ElevatorDirection.UP and destination <= person.home_floor:
                continue
            if current_direction == ElevatorDirection.DOWN and destination >= person.home_floor:
                continue
            # Separar solicitudes de trabajo
            if activity == ActivityType.GOING_TO_WORK:
                work_requests.append((person, destination, activity))
            elif not prioritize_work:
                compatible_requests.append((person, destination, activity))
        # Si hay solicitudes de trabajo, solo devolver esas
        return work_requests if work_requests else compatible_requests[:self.max_capacity - len(self.passengers)]

    def remove_request(self, person: Person):
        """Elimina una solicitud específica de la cola."""
        self.requests = [(p, pid, pr, d, a) for p, pid, pr, d, a in self.requests if pr != person]

    def add_passenger(self, person: Person, destination: int, activity: ActivityType) -> float:
        """Añade un pasajero al ascensor, simulando tiempo de entrada (0.5s)."""
        if len(self.passengers) < self.max_capacity:
            person.state = PersonState.IN_ELEVATOR
            person.destination_floor = destination
            person.activity = activity
            person.request_time = None  # Limpiar tiempo de solicitud
            self.passengers.append(person)
            time.sleep(0.5)
            return 0.5
        return 0.0

    def remove_passenger(self, person: Person) -> float:
        """Saca a un pasajero, simulando tiempo de salida (0.5s)."""
        if person in self.passengers:
            self.passengers.remove(person)
            person.state = PersonState.OUT
            # Cambiar a WORKING si llegó al piso 0 con GOING_TO_WORK
            if person.activity == ActivityType.GOING_TO_WORK and person.destination_floor == 0:
                person.activity = ActivityType.WORKING
            person.destination_floor = None
            person.request_time = None
            time.sleep(0.5)
            return 0.5
        return 0.0

    def move_to(self, target_floor: int, update_callback: Callable[['Elevator', List[str]], None], building: 'Building') -> float:
        """Mueve el ascensor al piso objetivo, recogiendo pasajeros en el camino."""
        total_time = 0.0
        new_events = []
        if self.current_floor != target_floor:
            new_events.append(f"Ascensor inicia movimiento desde piso {self.current_floor} a {target_floor}")
            time.sleep(2.0)
            total_time += 2.0
            step = 1 if target_floor > self.current_floor else -1
            self.direction = ElevatorDirection.UP if step > 0 else ElevatorDirection.DOWN
            prioritize_work = any(p.activity == ActivityType.GOING_TO_WORK for p in self.passengers)
            while self.current_floor != target_floor:
                self.current_floor += step
                new_events.append(f"Ascensor en piso {self.current_floor}")
                # Verificar si hay pasajeros para recoger en este piso
                compatible_requests = self.get_compatible_requests(self.direction, prioritize_work, building.floors)
                for person, destination, activity in compatible_requests:
                    if person.home_floor == self.current_floor and person in building.floors[self.current_floor]:
                        entry_time = self.add_passenger(person, destination, activity)
                        if entry_time > 0:
                            new_events.append(f"Persona {person.id} entra al ascensor en piso {self.current_floor} ({activity.value})")
                            building.floors[self.current_floor].remove(person)
                            self.remove_request(person)
                # Verificar si hay pasajeros para dejar en este piso
                for passenger in self.passengers[:]:
                    if passenger.destination_floor == self.current_floor:
                        exit_time = self.remove_passenger(passenger)
                        if exit_time > 0:
                            new_events.append(f"Persona {passenger.id} sale del ascensor en piso {self.current_floor} ({passenger.activity.value})")
                            building.floors[self.current_floor].append(passenger)
                update_callback(self, new_events)
                time.sleep(2.0)
                total_time += 2.0
            time.sleep(2.0)
            total_time += 2.0
        self.direction = ElevatorDirection.IDLE
        update_callback(self, new_events)
        return total_time

# Clase Display
class Display(ABC):
    """Interfaz para mostrar el estado del edificio y el ascensor."""
    @abstractmethod
    def update(self, building: 'Building', elevator: Elevator, events: List[str]) -> None:
        pass

# Clase CursesDisplay
class CursesDisplay(Display):
    """Muestra el estado del edificio y el ascensor en una interfaz tabular."""
    def __init__(self, stdscr):
        """Inicializa la ventana de curses."""
        self.stdscr = stdscr
        curses.curs_set(0)
        self.stdscr.timeout(100)

    def update(self, building: 'Building', elevator: Elevator, events: List[str]) -> None:
        """Actualiza la pantalla con el estado del edificio y ascensor en una tabla."""
        self.stdscr.clear()
        height, width = self.stdscr.getmaxyx()
        # Dividir pantalla en dos columnas
        col_width = width // 2
        self.stdscr.addstr(0, 0, "Pisos".ljust(col_width - 1) + "|" + "Ascensor".center(col_width - 1))
        self.stdscr.addstr(1, 0, "-" * (width - 1))
        for floor in range(building.total_floors - 1, -1, -1):
            waiting = [f"{p.id}({p.get_waiting_display()})" for p in building.floors[floor] if p.state == PersonState.WAITING]
            floor_str = f"Piso {floor}: {waiting}"[:col_width - 2]
            elevator_str = ""
            if elevator.current_floor == floor:
                passenger_info = [f"{p.id}({p.activity.value})" for p in elevator.passengers]
                elevator_str = f"{passenger_info} ({elevator.direction.value})"[:col_width - 2]
            self.stdscr.addstr(building.total_floors - floor + 2, 0, floor_str.ljust(col_width - 1) + "|" + elevator_str.rjust(col_width - 1))
        event_start = building.total_floors + 4
        self.stdscr.addstr(event_start, 0, "Eventos:")
        for i, event in enumerate(events[-5:]):
            if event_start + i + 1 < height:
                self.stdscr.addstr(event_start + i + 1, 0, event[:width - 1])
        queue_status = f"Cola: {[(p, d, a.value, prio) for prio, _, p, d, a in elevator.requests]}"
        self.stdscr.addstr(event_start + 6, 0, queue_status[:width - 1])
        self.stdscr.refresh()

# Clase Building
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
        # Inicializar personas solo en pisos 1 a total_floors-1
        for floor in range(1, total_floors):
            num_people = random.randint(1, 5)
            for _ in range(num_people):
                self.floors[floor].append(Person(id=self.person_id, home_floor=floor))
                self.person_id += 1

    def calculate_priority(self, person: Person, destination: int, activity: ActivityType) -> float:
        """Calcula la prioridad de una solicitud basada en distancia, actividad y tiempo de espera."""
        distance = abs(self.elevator.current_floor - person.home_floor)
        activity_weights = {
            ActivityType.GOING_TO_WORK: 1,
            ActivityType.WORKING: 4,  # No afecta prioridad, ya que no genera solicitudes
            ActivityType.SHOPPING: 2,
            ActivityType.PARTY: 3,
            ActivityType.EAT: 4,
            ActivityType.SLEEP: 4,
            ActivityType.ENTERTAIN: 4
        }
        activity_priority = activity_weights.get(activity, 4)
        # Reducir prioridad según tiempo de espera (menor valor = mayor prioridad)
        wait_time = (time.time() - person.request_time) if person.request_time else 0
        wait_factor = wait_time * 0.1  # Reducir 0.1 por segundo de espera
        return distance * 10 + activity_priority - wait_factor

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
                    prioritize_work = activity == ActivityType.GOING_TO_WORK or any(p.activity == ActivityType.GOING_TO_WORK for p in self.elevator.passengers)
                    new_events.append(f"Ascensor atiende a persona {person.id} en piso {person.home_floor} ({activity.value}, prioridad {self.calculate_priority(person, dest_floor, activity):.1f})")
                    # Mover al piso de la persona y procesar recogidas intermedias
                    self.elevator.move_to(person.home_floor, self.update_display, self)
                    entry_time = self.elevator.add_passenger(person, dest_floor, activity)
                    if entry_time > 0:
                        new_events.append(f"Persona {person.id} entra al ascensor en piso {person.home_floor} ({activity.value})")
                        if person in self.floors[person.home_floor]:  # Verificación adicional
                            self.floors[person.home_floor].remove(person)
                        self.elevator.remove_request(person)
                    else:
                        new_events.append(f"Ascensor lleno, persona {person.id} espera")
                        self.elevator.add_request(person, dest_floor, activity, self.calculate_priority(person, dest_floor, activity))
                    # Procesar pasajeros que se bajan en el destino
                    for passenger in self.elevator.passengers[:]:
                        if passenger.destination_floor == self.elevator.current_floor:
                            exit_time = self.elevator.remove_passenger(passenger)
                            if exit_time > 0:
                                new_events.append(f"Persona {passenger.id} sale del ascensor en piso {self.elevator.current_floor} ({passenger.activity.value})")
                                self.floors[self.elevator.current_floor].append(passenger)
                    # Mover al destino principal
                    self.elevator.move_to(dest_floor, self.update_display, self)
                    for passenger in self.elevator.passengers[:]:
                        if passenger.destination_floor == self.elevator.current_floor:
                            exit_time = self.elevator.remove_passenger(passenger)
                            if exit_time > 0:
                                new_events.append(f"Persona {passenger.id} sale del ascensor en piso {self.elevator.current_floor} ({passenger.activity.value})")
                                self.floors[self.elevator.current_floor].append(passenger)
                else:
                    new_events.append(f"Solicitud de persona {person.id} ignorada: estado inválido o no en piso {person.home_floor}")
                    self.elevator.remove_request(person)  # Limpiar solicitudes inválidas

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
