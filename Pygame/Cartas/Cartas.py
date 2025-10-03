from enum import Enum
from typing import List, Set, Optional, Dict, Union
import sys

# Tipos de cartas
class CardType(Enum):
    MATERIAL = "Material"
    DESIGN = "Design"
    FACTORY = "Factory"

# Clase base para todas las cartas
class Card:
    def __init__(self, card_type: CardType, card_id: int, name: str):
        self.card_type = card_type
        self.card_id = card_id
        self.name = name

    def __str__(self):
        return f"{self.card_type.value}: {self.name} (ID: {self.card_id})"

# Clase Material
class Material(Card):
    def __init__(self, card_id: int, name: str, tags: Set[str], composition: Set[int]):
        super().__init__(CardType.MATERIAL, card_id, name)
        self.tags = tags  # Etiquetas como 'metal', 'wood', etc.
        self.composition = composition  # IDs de materiales básicos que lo componen

    def consume_r(self) -> None:
        pass  # Marca el material como consumido

    def consume_m(self) -> None:
        pass  # Similar, sin objetos intermedios

    def has_component_r(self, component_id: int) -> bool:
        return component_id in self.composition

    def has_component_m(self, component_id: int) -> bool:
        for comp_id in self.composition:
            if comp_id == component_id:
                return True
        return False

# Clase Design
class Design(Card):
    def __init__(self, card_id: int, name: str, input_tags: Set[str], output_material_id: int, is_recycle: bool = False):
        super().__init__(CardType.DESIGN, card_id, name)
        self.input_tags = input_tags  # Etiquetas requeridas para insumos (e.g., {'wood', 'metal'})
        self.output_material_id = output_material_id
        self.is_recycle = is_recycle  # Indica si es diseño de reciclaje

    def is_compatible_r(self, factory: 'Factory') -> bool:
        return self.card_id in factory.allowed_design_ids

    def is_compatible_m(self, factory: 'Factory') -> bool:
        for allowed_id in factory.allowed_design_ids:
            if allowed_id == self.card_id:
                return True
        return False

    def check_materials_r(self, materials: List[Material]) -> bool:
        provided_tags = set()
        for mat in materials:
            provided_tags.update(mat.tags)
        return self.input_tags.issubset(provided_tags)

    def check_materials_m(self, materials: List[Material]) -> bool:
        required = list(self.input_tags)
        for mat in materials:
            for tag in mat.tags:
                if tag in required:
                    required.remove(tag)
                if not required:
                    return True
        return False

# Clase Factory
class Factory(Card):
    def __init__(self, card_id: int, name: str, allowed_design_ids: Set[int], is_recycling: bool = False):
        super().__init__(CardType.FACTORY, card_id, name)
        self.allowed_design_ids = allowed_design_ids
        self.is_recycling = is_recycling

    def produce_r(self, design: Design, materials: List[Material], material_db: Dict[int, tuple[str, Set[str], Set[int]]], consume: bool = True) -> Optional[Material]:
        if not design.is_compatible_r(self) or not design.check_materials_r(materials):
            return None
        if not self.is_recycling:
            if consume:
                for material in materials:
                    material.consume_r()
            name, tags, comp = material_db.get(design.output_material_id, ("Unknown", set(), set()))
            # Composición del nuevo material es la unión de las composiciones de insumos
            new_comp = set()
            for mat in materials:
                new_comp.update(mat.composition)
            return Material(design.output_material_id, name, tags, new_comp)
        else:
            if len(materials) != 1 or not design.is_recycle:
                return None
            input_material = materials[0]
            if input_material.has_component_r(design.output_material_id):
                if consume:
                    input_material.consume_r()
                name, tags, comp = material_db.get(design.output_material_id, ("Unknown", set(), set()))
                return Material(design.output_material_id, name, tags, comp)
        return None

    def produce_m(self, design: Design, materials: List[Material], material_db: Dict[int, tuple[str, Set[str], Set[int]]], consume: bool = True) -> Optional[Material]:
        if not design.is_compatible_m(self) or not design.check_materials_m(materials):
            return None
        if not self.is_recycling:
            if consume:
                for material in materials:
                    material.consume_m()
            name, tags, comp = material_db.get(design.output_material_id, ("Unknown", set(), set()))
            new_comp = set()
            for mat in materials:
                for c in mat.composition:
                    new_comp.add(c)
            return Material(design.output_material_id, name, tags, new_comp)
        else:
            if len(materials) != 1 or not design.is_recycle:
                return None
            input_material = materials[0]
            if input_material.has_component_m(design.output_material_id):
                if consume:
                    input_material.consume_m()
                name, tags, comp = material_db.get(design.output_material_id, ("Unknown", set(), set()))
                return Material(design.output_material_id, name, tags, comp)
        return None

    def recycle_r(self, designs: List[Design], material: Material, material_db: Dict[int, tuple[str, Set[str], Set[int]]]) -> List[Material]:
        if not self.is_recycling:
            return []
        result = []
        for design in designs:
            if design.is_compatible_r(self):
                produced = self.produce_r(design, [material], material_db, consume=False)
                if produced:
                    result.append(produced)
        if result:
            material.consume_r()
        return result

    def recycle_m(self, designs: List[Design], material: Material, material_db: Dict[int, tuple[str, Set[str], Set[int]]]) -> List[Material]:
        if not self.is_recycling:
            return []
        result = []
        for design in designs:
            if design.is_compatible_m(self):
                produced = self.produce_m(design, [material], material_db, consume=False)
                if produced:
                    result.append(produced)
        if result:
            material.consume_m()
        return result

# Clase Player
class Player:
    def __init__(self, name: str):
        self.name = name
        self.cards: List[Card] = []  # Lista de cartas del jugador

    def add_card(self, card: Card) -> None:
        self.cards.append(card)

    def remove_card(self, card: Card) -> None:
        self.cards.remove(card)

    def find_card_r(self, name: str) -> Optional[Card]:
        return next((c for c in self.cards if c.name.lower() == name.lower()), None)

    def find_card_m(self, name: str) -> Optional[Card]:
        for c in self.cards:
            if c.name.lower() == name.lower():
                return c
        return None

    def show_cards(self) -> str:
        return ", ".join(str(c) for c in self.cards)

# Función para parsear input (versión rápida)
def parse_input_r(input_str: str) -> List[str]:
    return [part.strip() for part in input_str.split('+')]

# Función para parsear input (versión memoria eficiente)
def parse_input_m(input_str: str) -> List[str]:
    parts = []
    current = ''
    for char in input_str:
        if char == '+':
            parts.append(current.strip())
            current = ''
        else:
            current += char
    if current:
        parts.append(current.strip())
    return parts

# Lógica del juego en modo texto
def text_game(player: Player, material_db: Dict[int, tuple[str, Set[str], Set[int]]], designs_db: Dict[str, Design], factories_db: Dict[str, Factory]):
    while True:
        print(f"Cartas de {player.name}: {player.show_cards()}")
        user_input = input("Ingresa acción (e.g., 'Madera + Acero + Diseño de Silla de Madera con Patas de Metal + Fabrica de muebles') o 'salir': ")
        if user_input.lower() == 'salir':
            break
        parts = parse_input_r(user_input)  # Usa _R por defecto
        if len(parts) < 3:
            print("Formato inválido. Usa 'Material(es) + Diseño(s) + Fábrica'")
            continue

        # Buscar todas las cartas mencionadas
        found_cards = []
        for part in parts:
            card = player.find_card_r(part)
            if card:
                found_cards.append(card)
            else:
                print(f"Carta '{part}' no encontrada.")
                break
        else:
            # Clasificar cartas
            materials = [c for c in found_cards if isinstance(c, Material)]
            designs = [c for c in found_cards if isinstance(c, Design)]
            factories = [c for c in found_cards if isinstance(c, Factory)]

            if len(factories) != 1 or len(designs) == 0 or len(materials) == 0:
                print("Debe haber exactamente una Fábrica, al menos un Diseño y al menos un Material.")
                continue

            factory = factories[0]

            if factory.is_recycling:
                if len(materials) != 1:
                    print("Reciclaje requiere exactamente un material.")
                    continue
                recycled = factory.recycle_r(designs, materials[0], material_db)
                if recycled:
                    player.remove_card(materials[0])
                    for rec in recycled:
                        player.add_card(rec)
                    print(f"Reciclado: {', '.join(str(r) for r in recycled)}")
                else:
                    print("Reciclaje fallido.")
            else:
                if len(designs) != 1:
                    print("Producción requiere exactamente un diseño.")
                    continue
                produced = factory.produce_r(designs[0], materials, material_db)
                if produced:
                    for mat in materials:
                        player.remove_card(mat)
                    player.add_card(produced)
                    print(f"Producido: {produced}")
                else:
                    print("Producción fallida.")

# Ejemplo de uso
def main():
    # Base de datos de materiales: {id: (name, tags, composition)}
    material_db = {
        1: ("Madera", {"wood"}, {1}),  # Básico
        2: ("Acero", {"metal"}, {2}),
        3: ("Silla de Madera", {"furniture"}, {1}),
        4: ("Silla de Madera con Patas de Metal", {"furniture"}, {1, 2})
    }

    # Diseños
    designs_db = {
        "Diseño de Silla de Madera con Patas de Metal": Design(1, "Diseño de Silla de Madera con Patas de Metal", {"wood", "metal"}, 4),
        "Diseño de Reciclaje de Madera": Design(2, "Diseño de Reciclaje de Madera", {"furniture"}, 1, is_recycle=True),
        "Diseño de Reciclaje de Metal": Design(3, "Diseño de Reciclaje de Metal", {"furniture"}, 2, is_recycle=True)
    }

    # Fábricas (ejemplo con una fábrica que no permite reciclaje de metal)
    factories_db = {
        "Fabrica de muebles": Factory(1, "Fabrica de muebles", {1}, is_recycling=False),
        "Fabrica de reciclaje": Factory(2, "Fabrica de reciclaje", {2, 3}, is_recycling=True),
        "Fabrica de reciclaje solo madera": Factory(3, "Fabrica de reciclaje solo madera", {2}, is_recycling=True)  # Solo permite reciclaje de madera
    }

    # Jugador
    player = Player("Jugador1")
    player.add_card(Material(1, "Madera", {"wood"}, {1}))
    player.add_card(Material(2, "Acero", {"metal"}, {2}))
    player.add_card(designs_db["Diseño de Silla de Madera con Patas de Metal"])
    player.add_card(factories_db["Fabrica de muebles"])
    player.add_card(designs_db["Diseño de Reciclaje de Madera"])
    player.add_card(designs_db["Diseño de Reciclaje de Metal"])
    player.add_card(factories_db["Fabrica de reciclaje"])
    player.add_card(factories_db["Fabrica de reciclaje solo madera"])

    text_game(player, material_db, designs_db, factories_db)

if __name__ == "__main__":
    main()
