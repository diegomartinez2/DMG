import random
import os
import sys

# Asegúrate de que ComfyUI esté en el path (ajusta si es necesario)
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Tablas predefinidas
personajes = [
    "un caballero valiente con armadura brillante", "un mago misterioso con un bastón brillante",
    "un hacker ciberpunk futurista", "un capitán pirata feroz",
    "una reina majestuosa con un vestido fluido", "un samurái solitario con katana",
    "un explorador en un dirigible steampunk", "un guerrero tribal con pintura de guerra"
]
escenas = [
    "de pie en un acantilado frente a un mar tormentoso", "en un mercado medieval bullicioso",
    "explorando una ciudad ciberpunk iluminada por neón", "dentro de un templo forestal antiguo",
    "en una nave espacial orbitando un planeta lejano", "en un desierto bajo un cielo estrellado",
    "vagando por un pantano brumoso", "en una cima montañosa nevada al amanecer",
    "en una selva vibrante con cascadas", "dentro de una cueva de cristal brillante"
]
estilos = [
    "realista", "anime", "pintura al óleo", "ciberpunk",
    "acuarela", "surrealista", "caricatura", "steampunk",
    "impresionista", "fotorrealista"
]
camaras = [
    "enfoque nítido, alta resolución", "enfoque suave, iluminación cinematográfica",
    "película granulada, aspecto vintage", "lente gran angular, perspectiva dramática",
    "profundidad de campo baja, bokeh", "4k ultra HD, hiperdetallado",
    "toma aérea con dron, vista expansiva", "lente macro, detalles intrincados"
]
fondos = [
    "una cadena montañosa brumosa", "un cielo al atardecer vibrante",
    "un paisaje urbano futurista con rascacielos", "un bosque exuberante con plantas brillantes",
    "una noche estrellada con auroras", "un castillo antiguo en ruinas",
    "un lago tranquilo rodeado de árboles otoñales", "un paisaje volcánico con lava fluyendo",
    "una playa serena con olas rompiendo", "un desierto azotado por el viento con dunas"
]
emociones = [
    "épico y heroico", "misterioso y espeluznante", "calmo y sereno",
    "intenso y dramático", "alegre y vibrante", "melancólico y sombrío",
    "etéreo y onírico", "caótico y turbulento"
]
aspect_ratios = ["16:9", "9:16", "1:1", "4:3", "3:4"]
versions = ["5", "1", "4"]
qualities = ["8", "16", "full"]

# Prompt negativo base
prompt_negativo_base = (
    "borroso, baja calidad, distorsionado, extremidades extras, extremidades faltantes, mala anatomía, "
    "iluminación pobre, sobresaturado, texto, marca de agua, baja resolución, pixelado, "
    "artefactos, dedos extras, cara deformada, proporciones antinaturales"
)

class GeneradorDePrompt:
    @classmethod
    def INPUT_TYPES(cls):
        return {
            "required": {
                "personaje": (["NINGUNO"] + personajes, {"default": "NINGUNO"}),
                "escena": (["NINGUNO"] + escenas, {"default": "NINGUNO"}),
                "estilo": (["NINGUNO"] + estilos, {"default": "NINGUNO"}),
                "camara": (["NINGUNO"] + camaras, {"default": "NINGUNO"}),
                "fondo": (["NINGUNO"] + fondos, {"default": "NINGUNO"}),
                "emocion": (["NINGUNO"] + emociones, {"default": "NINGUNO"}),
                "incluir_personaje": ("BOOLEAN", {"default": True}),
                "ar": (["NINGUNO"] + aspect_ratios, {"default": "NINGUNO"}),
                "v": (["NINGUNO"] + versions, {"default": "NINGUNO"}),
                "q": (["NINGUNO"] + qualities, {"default": "NINGUNO"}),
                "usar_aleatorio": ("BOOLEAN", {"default": True}),
                "seed": ("INT", {"default": 0, "min": 0, "max": 0xffffffffffffffff}),
            }
        }

    RETURN_TYPES = ("STRING", "STRING", "STRING", "STRING", "STRING")
    RETURN_NAMES = ("prompt_positivo", "prompt_negativo", "relacion_de_aspecto", "version_sampler", "calidad")
    FUNCTION = "generar"
    CATEGORY = "Utilidades/Prompts"

    def generar(self, personaje, escena, estilo, camara, fondo, emocion, incluir_personaje, ar, v, q, usar_aleatorio, seed):
        random.seed(seed)  # Usa seed para reproducibilidad

        # Función helper para obtener valor o aleatorio
        def get_or_random(val, options):
            if val == "NINGUNO" and usar_aleatorio:
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
            prompt_positivo = (
                f"{personaje}, {escena}, estilo {estilo}, {camara}, "
                f"fondo con {fondo}, atmósfera {emocion}, altamente detallado, colores vibrantes"
            )
        else:
            prompt_positivo = (
                f"un paisaje detallado de {escena}, estilo {estilo}, {camara}, "
                f"fondo con {fondo}, atmósfera {emocion}, altamente detallado, colores vibrantes, "
                f"detalles ambientales ricos, iluminación dinámica"
            )

        prompt_negativo = prompt_negativo_base

        # Formatear ar como "width:height"
        ar_formatted = self.format_ar(ar_val)

        return (prompt_positivo, prompt_negativo, ar_formatted, v_val, q_val)

    def format_ar(self, ar_str):
        ratios = {"16:9": "1024:576", "9:16": "576:1024", "1:1": "512:512", "4:3": "768:576", "3:4": "576:768"}
        return ratios.get(ar_str, "512:512")

# Registro del nodo
NODE_CLASS_MAPPINGS = {
    "GeneradorDePrompt": GeneradorDePrompt
}

NODE_DISPLAY_NAME_MAPPINGS = {
    "GeneradorDePrompt": "Generador de Prompt"
}
