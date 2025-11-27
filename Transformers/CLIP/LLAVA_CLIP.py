# pip install llama-index-llms-ollama llama-index-embeddings-huggingface pillow
from llama_index.multi_modal_llms.ollama import OllamaMultiModal
from PIL import Image

# Necesitas tener Ollama instalado y el modelo descargado:
# 1. Instala Ollama: https://ollama.com
# 2. En terminal: ollama pull llava:13b   (o llava:34b si tienes mucha VRAM)

mm_llm = OllamaMultiModal(model="llava:13b", temperature=0.1)

# Cargar tu imagen
imagen = Image.open("tu_imagen.jpg")   # o "feynman.png", "espectro.png", etc.

# Preguntar lo que quieras
respuesta = mm_llm.complete(
    prompt="Describe esta imagen con todo detalle científico posible. Si es un gráfico o diagrama, extrae datos y explica qué representa.",
    image_documents=[imagen],
)

print(respuesta)
