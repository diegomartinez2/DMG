import random

# Tablas predefinidas para entradas aleatorias
personajes = [
    "a brave knight in shining armor", "a mysterious wizard with a glowing staff",
    "a futuristic cyberpunk hacker", "a fierce pirate captain",
    "a regal queen in a flowing gown", "a lone samurai with a katana"
]
escenas = [
    "standing on a cliff overlooking a stormy sea", "in a bustling medieval marketplace",
    "exploring a neon-lit cyberpunk city", "inside an ancient forest temple",
    "on a spaceship orbiting a distant planet", "in a desert under a starry sky"
]
estilos = [
    "realistic", "anime", "oil painting", "cyberpunk",
    "watercolor", "surrealist", "cartoon", "steampunk"
]
camaras = [
    "sharp focus, high resolution", "soft focus, cinematic lighting",
    "grainy film, vintage look", "wide-angle lens, dramatic perspective",
    "shallow depth of field, bokeh", "4k ultra HD, hyper-detailed"
]
fondos = [
    "a misty mountain range", "a vibrant sunset sky",
    "a futuristic cityscape with skyscrapers", "a lush forest with glowing plants",
    "a starry night with auroras", "a ruined ancient castle"
]
emociones = [
    "epic and heroic", "mysterious and eerie", "calm and serene",
    "intense and dramatic", "joyful and vibrant", "melancholic and somber"
]

# Prompt negativo estándar para Stable Diffusion
prompt_negativo_base = (
    "blurry, low quality, distorted, extra limbs, missing limbs, bad anatomy, "
    "poor lighting, oversaturated, text, watermark, low resolution, pixelated, "
    "artifacts, extra fingers, deformed face, unnatural proportions"
)

def obtener_entrada_usuario(mensaje, opciones, campo):
    entrada = input(mensaje).strip()
    if not entrada:
        print(f"No se proporcionó {campo}, seleccionando uno aleatoriamente...")
        return random.choice(opciones)
    return entrada

def generar_prompt():
    print("Generador de Prompts para Stable Diffusion 1.5")
    print("Responde las siguientes preguntas. Presiona Enter para seleccionar un valor aleatorio.\n")

    # Solicitar entradas al usuario
    personaje = obtener_entrada_usuario(
        "Describe al personaje (ej. 'un caballero medieval', 'una astronauta'): ",
        personajes, "personaje"
    )
    escena = obtener_entrada_usuario(
        "Describe la escena o lugar (ej. 'en un bosque encantado', 'en una ciudad futurista'): ",
        escenas, "escena"
    )
    estilo = obtener_entrada_usuario(
        "Especifica el estilo artístico (ej. 'realista', 'anime', 'pintura al óleo'): ",
        estilos, "estilo artístico"
    )
    camara = obtener_entrada_usuario(
        "Describe el tipo de cámara o efecto visual (ej. 'enfoque nítido', 'granulado'): ",
        camaras, "tipo de cámara"
    )
    fondo = obtener_entrada_usuario(
        "Describe el fondo (ej. 'montañas brumosas', 'cielo estrellado'): ",
        fondos, "fondo"
    )
    emocion = obtener_entrada_usuario(
        "Describe la emoción o atmósfera (ej. 'épico', 'sombrío', 'alegre'): ",
        emociones, "emoción"
    )

    # Construir el prompt positivo
    prompt_positivo = (
        f"{personaje}, {escena}, {estilo} style, {camara}, "
        f"background with {fondo}, {emocion} atmosphere, highly detailed, vibrant colors"
    )

    # Mostrar resultados
    print("\n=== Prompt Generado ===")
    print("Prompt Positivo:", prompt_positivo)
    print("Prompt Negativo:", prompt_negativo_base)

    return prompt_positivo, prompt_negativo_base

def main():
    while True:
        generar_prompt()
        continuar = input("\n¿Quieres generar otro prompt? (s/n): ").strip().lower()
        if continuar != 's':
            print("¡Gracias por usar el generador de prompts!")
            break

if __name__ == "__main__":
    main()
