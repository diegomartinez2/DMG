import random

def generate_sci_fi_names(count=20):
    # Apellidos más comunes globalmente (India, China, LatAm, Anglo, etc.)
    last_names = [
        "Singh", "Wang", "Devi", "Li", "Zhang", "García", "Khan", "Chen",
        "Rodríguez", "Smith", "Liu", "Kumar", "Nguyen", "Hernández", "Müller",
        "Suzuki", "Chatterjee", "Jones", "Ahmed", "Kim", "Sato", "Ivanov"
    ]

    # Nombres por "Modas" o Eras
    eras = {
        "Legado (Tradicional 2020)": [
            "Mateo", "Sofía", "Liam", "Emma", "Noah", "Olivia", "Lucas", "Isabella",
            "Santiago", "Mia", "Ethan", "Aria", "Sebastian", "Zoe"
        ],
        "Tecno-Optimista (Conceptuales)": [
            "Nova", "Zenit", "Átomo", "Caelum", "Data", "Índigo", "Vortex", "Soma",
            "Aura", "Éter", "Órbita", "Silicio", "Vector", "Helix"
        ],
        "Amalgama (Sintéticos/Breves)": [
            "Ixo", "Kly", "Veda", "Runa", "Jyx", "Mio", "Oru", "Tal",
            "Nyx", "Zel", "Aks", "Eon", "Vyr", "Luz"
        ],
        "Neo-Clásicos (Latín/Griego recuperado)": [
            "Aethel", "Cora", "Dante", "Gaya", "Helios", "Io", "Juno", "Lira",
            "Merope", "Nox", "Orion", "Pax", "Talia", "Xeno"
        ]
    }

    generated_list = []

    for _ in range(count):
        era_name = random.choice(list(eras.keys()))
        first_name = random.choice(eras[era_name])
        last_name = random.choice(last_names)

        # Ocasionalmente un nombre doble o compuesto para añadir variedad
        if random.random() > 0.85:
            first_name += "-" + random.choice(eras[era_name])

        generated_list.append({
            "Nombre": f"{first_name} {last_name}",
            "Era": era_name
        })

    return generated_list

# Ejecución del generador
if __name__ == "__main__":
    nombres = generate_sci_fi_names(25)

    print(f"{'PERSONAJE':<30} | {'ERA/MODA DE NACIMIENTO'}")
    print("-" * 60)
    for p in nombres:
        print(f"{p['Nombre']:<30} | {p['Era']}")
