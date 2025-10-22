import random
import os
import sys
import importlib.util

# Asegúrate de que ComfyUI esté en el path (ajusta si es necesario)
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from comfy import model_management  # Para manejo de modelos si se necesita

# Tablas predefinidas (mismas que en el script original)
personajes = [
    "a brave knight in shining armor", "a mysterious wizard with a glowing staff",
    "a futuristic cyberpunk hacker", "a fierce pirate captain",
    "a regal queen in a flowing gown", "a lone samurai with a katana",
    "an explorer in a steampunk airship", "a tribal warrior with war paint"
]
escenas = [
    "standing on a cliff overlooking a stormy sea", "in a bustling medieval marketplace",
    "exploring a neon-lit cyberpunk city", "inside an ancient forest temple",
    "on a spaceship orbiting a distant planet", "in a desert under a starry sky",
    "wandering through a misty swamp", "on a snowy mountain peak at dawn",
    "in a vibrant jungle with cascading waterfalls", "inside a glowing crystal cave"
]
estilos = [
    "realistic", "anime", "oil painting", "cyberpunk",
    "watercolor", "surrealist", "cartoon", "steampunk",
    "impressionist", "photorealistic"
]
camaras = [
    "sharp focus, high resolution", "soft focus, cinematic lighting",
    "grainy film, vintage look", "wide-angle lens, dramatic perspective",
    "shallow depth of field, bokeh", "4k ultra HD, hyper-detailed",
    "aerial drone shot, expansive view", "macro lens, intricate details"
]
fondos = [
    "a misty mountain range", "a vibrant sunset sky",
    "a futuristic cityscape with skyscrapers", "a lush forest with glowing plants",
    "a starry night with auroras", "a ruined ancient castle",
    "a tranquil lake surrounded by autumn trees", "a volcanic landscape with flowing lava",
    "a serene beach with crashing waves", "a windswept desert with sand dunes"
]
emociones = [
    "epic and heroic", "mysterious and eerie", "calm and serene",
    "intense and dramatic", "joyful and vibrant", "melancholic and somber",
    "ethereal and dreamlike", "chaotic and turbulent"
]
aspect_ratios = ["16:9", "9:16", "1:1", "4:3", "3:4"]
versions = ["5", "1", "4"]
qualities = ["8", "16", "full"]

# Prompt negativo base
prompt_negativo_base = (
    "blurry, low quality, distorted, extra limbs, missing limbs, bad anatomy, "
    "poor lighting, oversaturated, text, watermark, low resolution, pixelated, "
    "artifacts, extra fingers, deformed face, unnatural proportions"
)

class PromptGeneratorNode:
    @classmethod
    def INPUT_TYPES(cls):
        return {
            "required": {
                "personaje": (["NONE"] + personajes, {"default": "NONE"}),
                "escena": (["NONE"] + escenas, {"default": "NONE"}),
                "estilo": (["NONE"] + estilos, {"default": "NONE"}),
                "camara": (["NONE"] + camaras, {"default": "NONE"}),
                "fondo": (["NONE"] + fondos, {"default": "NONE"}),
                "emocion": (["NONE"] + emociones, {"default": "NONE"}),
                "incluir_personaje": ("BOOLEAN", {"default": True}),
                "ar": (["NONE"] + aspect_ratios, {"default": "NONE"}),
                "v": (["NONE"] + versions, {"default": "NONE"}),
                "q": (["NONE"] + qualities, {"default": "NONE"}),
                "use_random": ("BOOLEAN", {"default": True}),  # Si True, usa random si NONE
                "seed": ("INT", {"default": 0, "min": 0, "max": 0xffffffffffffffff}),
            }
        }

    RETURN_TYPES = ("STRING", "STRING", "STRING", "STRING", "STRING")  # positive, negative, ar_formatted, v, q
    RETURN_NAMES = ("positive_prompt", "negative_prompt", "aspect_ratio", "sampler_version", "quality")
    FUNCTION = "generate"
    CATEGORY = "PromptUtils"  # Categoría en el menú de nodos

    def generate(self, personaje, escena, estilo, camara, fondo, emocion, incluir_personaje, ar, v, q, use_random, seed):
        random.seed(seed)  # Usa seed para reproducibilidad

        # Función helper para obtener valor o random
        def get_or_random(val, options):
            if val == "NONE" or (use_random and val == "NONE"):
                return random.choice(options)
            return val

        # Obtener valores
        personaje = get_or_random(personaje, personajes) if incluir_personaje else ""
        escena = get_or_random(escena, escenas)
        estilo = get_or_random(estilo, estilos)
        camara = get_or_random(camara, camaras)
        fondo = get_or_random(fondo, fondos)
        emocion = get_or_random(emocion, emociones)
        ar_val = get_or_random(ar, aspect_ratios)
        v_val = get_or_random(v, versions)
        q_val = get_or_random(q, qualities)

        # Construir prompt positivo
        if incluir_personaje:
            positive = (
                f"{personaje}, {escena}, {estilo} style, {camara}, "
                f"background with {fondo}, {emocion} atmosphere, highly detailed, vibrant colors"
            )
        else:
            positive = (
                f"a detailed landscape of {escena}, {estilo} style, {camara}, "
                f"background with {fondo}, {emocion} atmosphere, highly detailed, vibrant colors, "
                f"rich environmental details, dynamic lighting"
            )

        negative = prompt_negativo_base

        # Formatear ar como "width:height" (ej. para nodo Empty Latent Image)
        ar_formatted = self.format_ar(ar_val)

        return (positive, negative, ar_formatted, v_val, q_val)

    def format_ar(self, ar_str):
        ratios = {"16:9": "1024:576", "9:16": "576:1024", "1:1": "512:512", "4:3": "768:576", "3:4": "576:768"}
        return ratios.get(ar_str, "512:512")

# Registro del nodo (se llama automáticamente al cargar)
NODE_CLASS_MAPPINGS = {
    "PromptGeneratorNode": PromptGeneratorNode
}

NODE_DISPLAY_NAME_MAPPINGS = {
    "PromptGeneratorNode": "Prompt Generator Node"
}
