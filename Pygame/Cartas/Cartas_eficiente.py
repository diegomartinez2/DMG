from enum import Enum
from typing import List, Optional, Dict
from collections import Counter
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
    def __init__(self, card_id: int, name: str, tags: frozenset[str], composition: List[int]):
        super().__init__(CardType.MATERIAL, card_id, name)
        self.tags = tags
        self.composition = Counter(composition)  # Contador para componentes

    def consume_r(self) -> None:
        pass

    def consume_m(self) -> None:
        pass

    def has_component_r(self, component_id: int) -> bool:
        return self.composition[component_id] > 0

    def has_component_m(self, component_id: int) -> bool:
        return self.composition.get(component_id, 0) > 0

    def add_components_r(self, components: List[int]) -> None:
        self.composition.update(components)

    def add_components_m(self, components: List[int]) -> None:
        for comp in components:
            self.composition[comp] += 1

    def remove_component_r(self, component_id: int) -> bool:
        if self.composition[component_id] > 0:
            self.composition[component_id] -= 1
            if self.composition[component_id] == 0:
                del self.composition[component_id]
            return True
        return False

    def remove_component_m(self, component_id: int) -> bool:
        count = self.composition.get(component_id, 0)
        if count > 0:
            self.composition[component_id] = count - 1
            if self.composition[component_id] == 0:
                del self.composition[component_id]
            return True
        return False

# Clase Design
class Design(Card):
    def __init__(self, card_id: int, name: str, input_tags: frozenset[str], output_material_id: int, is_recycle: bool = False):
        super().__init__(CardType.DESIGN, card_id, name)
        self.input_tags = input_tags
        self.output_material_id = output_material_id
        self.is_recycle = is_recycle

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
    def __init__(self, card_id: int, name: str, allowed_design_ids: frozenset[int], is_recycling: bool = False):
        super().__init__(CardType.FACTORY, card_id, name)
        self.allowed_design_ids = allowed_design_ids
        self.is_recycling = is_recycling

    def produce_r(self, design: Design, materials: List[Material], material_db: Dict[int, tuple[str, frozenset[str], List[int]]], consume: bool = True) -> Optional[Material]:
        if not design.is_compatible_r(self) or not design.check_materials_r(materials):
            return None
        if not self.is_recycling:
            if consume:
                for material in materials:
                    material.consume_r()
            name, tags, _ = material_db.get(design.output_material_id, ("Unknown", frozenset(), []))
            new_comp = []
            for mat in materials:
                new_comp.extend(mat.composition.elements())
            return Material(design.output_material_id, name, tags, new_comp)
        else:
            if len(materials) != 1 or not design.is_recycle:
                return None
            input_material = materials[0]
            if input_material.has_component_r(design.output_material_id):
                if consume:
                    input_material.consume_r()
                name, tags, comp = material_db.get(design.output_material_id, ("Unknown", frozenset(), []))
                return Material(design.output_material_id, name, tags, comp)
        return None

    def produce_m(self, design: Design, materials: List[Material], material_db: Dict[int, tuple[str, frozenset[str], List[int]]], consume: bool = True) -> Optional[Material]:
        if not design.is_compatible_m(self) or not design.check_materials_m(materials):
            return None
        if not self.is_recycling:
            if consume:
                for material in materials:
                    material.consume_m()
            name, tags, _ = material_db.get(design.output_material_id, ("Unknown", frozenset(), []))
            new_comp = []
            for mat in materials:
                for comp_id, count in mat.composition.items():
                    new_comp.extend([comp_id] * count)
            return Material(design.output_material_id, name, tags, new_comp)
        else:
            if len(materials) != 1 or not design.is_recycle:
                return None
            input_material = materials[0]
            if input_material.has_component_m(design.output_material_id):
                if consume:
                    input_material.consume_m()
                name, tags, comp = material_db.get(design.output_material_id, ("Unknown", frozenset(), []))
                return Material(design.output_material_id, name, tags, comp)
        return None

    def recycle_r(self, designs: List[Design], material: Material, material_db: Dict[int, tuple[str, frozenset[str], List[int]]]) -> List[Material]:
        if not self.is_recycling:
            return []
        result = []
        leftover = Counter(material.composition)
        for design in designs:
            if design.is_compatible_r(self) and material.has_component_r(design.output_material_id):
                produced = self.produce_r(design, [material], material_db, consume=False)
                if produced and material.remove_component_r(design.output_material_id):
                    result.append(produced)
                    leftover[design.output_material_id] -= 1
                    if leftover[design.output_material_id] == 0:
                        del leftover[design.output_material_id]
        if leftover:
            result.append(Material(999, "Basura", frozenset(["waste"]), list(leftover.elements())))
        if result:
            material.consume_r()
        return result

    def recycle_m(self, designs: List[Design], material: Material, material_db: Dict[int, tuple[str, frozenset[str], List[int]]]) -> List[Material]:
        if not self.is_recycling:
            return []
        result = []
        leftover = Counter()
        for comp_id, count in material.composition.items():
            leftover[comp_id] = count
        for design in designs:
            if design.is_compatible_m(self) and material.has_component_m(design.output_material_id):
                produced = self.produce_m(design, [material], material_db, consume=False)
                if produced and material.remove_component_m(design.output_material_id):
                    result.append(produced)
                    leftover[design.output_material_id] -= 1
                    if leftover[design.output_material_id] == 0:
                        del leftover[design.output_material_id]
        if leftover:
            result.append(Material(999, "Basura", frozenset(["waste"]), list(leftover.elements())))
        if result:
            material.consume_m()
        return result

# Clase Player
class Player:
    def __init__(self, name: str):
        self.name = name
        self.cards: List[Card] = []

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

    def add_basura_r(self, new_basura: Material) -> None:
        existing_basura = self.find_card_r("Basura")
        if existing_basura:
            existing_basura.add_components_r(new_basura.composition.elements())
        else:
            self.add_card(new_basura)

    def add_basura_m(self, new_basura: Material) -> None:
        existing_basura = self.find_card_m("Basura")
        if existing_basura:
            existing_basura.add_components_m(new_basura.composition.elements())
        else:
            self.add_card(new_basura)

# Función para parsear input (versión rápida)
def parse_input_r(input_str: str) -> List[str]:
    return [part.strip() for part in input_str.split('+')]

# Función para parsear input (versión memoria eficiente)
def parse_input_m(input_str: str) -> List[str]:
    parts = []
    current = ''
    for char in input_str:
        if char == '+':
            if current.strip():
                parts.append(current.strip())
            current = ''
        else:
            current += char
    if current.strip():
        parts.append(current.strip())
    return parts

# Lógica del juego en modo texto
def text_game(player: Player, material_db: Dict[int, tuple[str, frozenset[str], List[int]]], designs_db: Dict[str, Design], factories_db: Dict[str, Factory]):
    while True:
        print(f"Cartas de {player.name}: {player.show_cards()}")
        user_input = input("Ingresa acción (e.g., 'Madera + Diseño de Silla de Madera con Patas de Metal + Acero + Fabrica de muebles') o 'salir': ")
        if user_input.lower() == 'salir':
            break
        parts = parse_input_r(user_input)
        if len(parts) < 3:
            print("Formato inválido. Usa al menos un Material, un Diseño y una Fábrica.")
            continue

        materials = []
        designs = []
        factories = []
        material_names = {name for id, (name, _, _) in material_db.items()}
        for part in parts:
            card = player.find_card_r(part)
            if not card:
                print(f"Carta '{part}' no encontrada.")
                break
            if part in designs_db:
                designs.append(card)
            elif part in factories_db:
                factories.append(card)
            elif part in material_names or part.lower().startswith("basura"):
                materials.append(card)
            else:
                print(f"Carta '{part}' no reconocida como Material, Diseño o Fábrica.")
                break
        else:
            if len(factories) != 1:
                print("Debe haber exactamente una Fábrica.")
                continue
            if len(designs) == 0:
                print("Debe haber al menos un Diseño.")
                continue
            if len(materials) == 0:
                print("Debe haber al menos un Material.")
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
                        if rec.card_id == 999:
                            player.add_basura_r(rec)
                        else:
                            player.add_card(rec)
                    print(f"Reciclado: {', '.join(str(r) for r in recycled)}")
                else:
                    print("Reciclaje fallido: ningún diseño compatible o material no válido.")
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
                    print("Producción fallida: diseño no compatible o materiales insuficientes.")

# Ejemplo de uso
def main():
    # Base de datos de materiales: {id: (name, tags, composition)}
    material_db = {
        1: ("Madera", frozenset(["wood"]), [1]),
        2: ("Acero", frozenset(["metal"]), [2]),
        3: ("Silla de Madera", frozenset(["furniture"]), [1]),
        4: ("Silla de Madera con Patas de Metal", frozenset(["furniture"]), [1, 2]),
        999: ("Basura", frozenset(["waste"]), [])
    }

    # Diseños
    designs_db = {
        "Diseño de Silla de Madera con Patas de Metal": Design(1, "Diseño de Silla de Madera con Patas de Metal", frozenset(["wood", "metal"]), 4),
        "Diseño de Reciclaje de Madera": Design(2, "Diseño de Reciclaje de Madera", frozenset(["furniture", "waste"]), 1, is_recycle=True),
        "Diseño de Reciclaje de Metal": Design(3, "Diseño de Reciclaje de Metal", frozenset(["furniture", "waste"]), 2, is_recycle=True)
    }

    # Fábricas
    factories_db = {
        "Fabrica de muebles": Factory(1, "Fabrica de muebles", frozenset([1]), is_recycling=False),
        "Fabrica de reciclaje": Factory(2, "Fabrica de reciclaje", frozenset([2, 3]), is_recycling=True),
        "Fabrica de reciclaje solo madera": Factory(3, "Fabrica de reciclaje solo madera", frozenset([2]), is_recycling=True)
    }

    # Jugador
    player = Player("Jugador1")
    player.add_card(Material(1, "Madera", frozenset(["wood"]), [1]))
    player.add_card(Material(2, "Acero", frozenset(["metal"]), [2]))
    player.add_card(Material(3, "Silla de Madera", frozenset(["furniture"]), [1]))
    player.add_card(Material(4, "Silla de Madera con Patas de Metal", frozenset(["furniture"]), [1, 2]))
    player.add_card(designs_db["Diseño de Silla de Madera con Patas de Metal"])
    player.add_card(designs_db["Diseño de Reciclaje de Madera"])
    player.add_card(designs_db["Diseño de Reciclaje de Metal"])
    player.add_card(factories_db["Fabrica de muebles"])
    player.add_card(factories_db["Fabrica de reciclaje"])
    player.add_card(factories_db["Fabrica de reciclaje solo madera"])

    # Simulación: combinar basuras
    print("Simulación: reciclar Silla de Madera con Patas de Metal solo para Madera")
    recycled = factories_db["Fabrica de reciclaje solo madera"].recycle_r(
        [designs_db["Diseño de Reciclaje de Madera"]],
        Material(4, "Silla de Madera con Patas de Metal", frozenset(["furniture"]), [1, 2]),
        material_db
    )
    for rec in recycled:
        player.add_basura_r(rec) if rec.card_id == 999 else player.add_card(rec)
    print(f"Reciclado: {', '.join(str(r) for r in recycled)}")
    print(f"Cartas de {player.name}: {player.show_cards()}")

    print("Simulación: reciclar otra Silla para añadir más Madera a Basura")
    recycled = factories_db["Fabrica de reciclaje solo madera"].recycle_r(
        [designs_db["Diseño de Reciclaje de Madera"]],
        Material(3, "Silla de Madera", frozenset(["furniture"]), [1]),
        material_db
    )
    for rec in recycled:
        player.add_basura_r(rec) if rec.card_id == 999 else player.add_card(rec)
    print(f"Reciclado: {', '.join(str(r) for r in recycled)}")
    print(f"Cartas de {player.name}: {player.show_cards()}")

    text_game(player, material_db, designs_db, factories_db)

if __name__ == "__main__":
    main()
