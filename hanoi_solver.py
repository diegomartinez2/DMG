class HanoiSolver:
    """
    Resuelve el puzle de las Torres de Hanoi para un número dado de discos.

    Esta clase aplica el Principio de Responsabilidad Única (SRP) al
    enfocarse únicamente en la lógica de resolución y la gestión del estado,
    separándola de cualquier lógica de visualización.

    La lógica de resolución implementa el principio KISS usando un algoritmo
    recursivo simple y elegante.
    """

    def __init__(self, num_disks: int):
        """
        Inicializa el solucionador.

        Args:
            num_disks (int): El número de discos a mover.
        """
        if num_disks < 1:
            raise ValueError("El número de discos debe ser al menos 1.")

        self.num_disks = num_disks
        self.moves = []  # Almacenará la lista de movimientos (tuplas)

        # El estado de las varillas (SRP: la clase gestiona su estado)
        # Usamos listas como pilas. Los discos se representan por números (tamaño)
        # El disco más pequeño (1) está al final de la lista 'A'.
        self.rods = {
            'A': list(range(num_disks, 0, -1)),  # Fuente
            'B': [],                             # Auxiliar
            'C': []                              # Destino
        }

    def solve(self) -> list[tuple[str, str]]:
        """
        Punto de entrada público para resolver el puzle.

        Aplica el OCP: la lógica interna está encapsulada, el usuario
        solo llama a solve().

        Returns:
            list[tuple[str, str]]: Una lista de tuplas, donde cada tupla
                                   representa un movimiento (desde_varilla, hasta_varilla).
        """
        # Limpiamos los movimientos por si se llama a solve() múltiples veces
        self.moves = []

        # Inicializamos el estado por si acaso
        self.rods = {
            'A': list(range(self.num_disks, 0, -1)),
            'B': [],
            'C': []
        }

        print(f"Estado inicial para {self.num_disks} discos:")
        self.print_rods()

        # Llamamos al ayudante recursivo privado
        self._move_disks_recursive(self.num_disks, 'A', 'C', 'B')

        print("\n¡Puzle resuelto!")
        self.print_rods()

        return self.moves

    def _move_disks_recursive(self, n: int, source: str, destination: str, auxiliary: str):
        """
        El núcleo del algoritmo recursivo (KISS).

        Mueve 'n' discos desde 'source' a 'destination' usando 'auxiliary'.

        Args:
            n (int): El número de discos a mover en este paso.
            source (str): La varilla de origen ('A', 'B', o 'C').
            destination (str): La varilla de destino ('A', 'B', o 'C').
            auxiliary (str): La varilla auxiliar ('A', 'B', o 'C').
        """
        # Caso base (KISS): Mover 1 disco es trivial
        if n == 1:
            self._move_single_disk(source, destination)
            return

        # Paso 1: Mover n-1 discos de Fuente a Auxiliar
        self._move_disks_recursive(n - 1, source, auxiliary, destination)

        # Paso 2: Mover el disco n (el más grande) de Fuente a Destino
        self._move_single_disk(source, destination)

        # Paso 3: Mover n-1 discos de Auxiliar a Destino
        self._move_disks_recursive(n - 1, auxiliary, destination, source)

    def _move_single_disk(self, source: str, destination: str):
        """
        Mueve físicamente un disco y registra el movimiento.

        Esta es una función de ayuda interna (SRP) que maneja la
        manipulación del estado y el registro de movimientos.
        """
        try:
            # 1. Quitar el disco de la parte superior de la varilla 'source'
            disk = self.rods[source].pop()

            # 2. Validar el movimiento (opcional pero buena práctica)
            if self.rods[destination] and self.rods[destination][-1] < disk:
                # Esto nunca debería ocurrir si la lógica es correcta,
                # pero es una buena salvaguarda.
                raise Exception(f"Movimiento inválido: Disco {disk} no puede ir sobre {self.rods[destination][-1]}")

            # 3. Poner el disco en la varilla 'destination'
            self.rods[destination].append(disk)

            # 4. Registrar el movimiento
            self.moves.append((source, destination))

            # print(f"Mover disco {disk} de {source} a {destination}")
            # self.print_rods() # Descomentar para depuración paso a paso

        except IndexError:
            print(f"Error: Se intentó mover desde una varilla vacía '{source}'")

    def print_rods(self):
        """Método de utilidad para imprimir el estado actual de las varillas."""
        print("---")
        for name, disks in self.rods.items():
            # Imprime los discos con el más grande (base) a la izquierda
            print(f"{name}: {disks}")
        print("---")

    def get_moves(self) -> list[tuple[str, str]]:
        """Retorna la lista de movimientos calculados."""
        return self.moves

# --- Ejemplo de uso ---
if __name__ == "__main__":

    NUM_DISCOS = 3

    # 1. Instanciamos nuestro módulo solucionador
    juego_hanoi = HanoiSolver(NUM_DISCOS)

    # 2. Llamamos al método público para resolver
    # La clase se encarga de toda la lógica interna.
    lista_de_movimientos = juego_hanoi.solve()

    # 3. Usamos los resultados
    print(f"\nTotal de movimientos: {len(lista_de_movimientos)}")

    # El módulo 'juego_hanoi' puede ser reutilizado. Por ejemplo,
    # una GUI podría llamar a get_moves() y animar cada paso.
    print("Secuencia de movimientos (Desde, Hasta):")
    for move in lista_de_movimientos:
        print(f"  {move[0]} -> {move[1]}")

    print("\n--- Probando con 4 discos ---")
    juego_4_discos = HanoiSolver(4)
    movimientos_4_discos = juego_4_discos.solve()
    print(f"\nTotal de movimientos para 4 discos: {len(movimientos_4_discos)}")
