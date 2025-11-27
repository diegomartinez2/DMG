# pip install clip torch torchvision pillow
import clip
import torch
from PIL import Image
import requests
from io import BytesIO

# Cargar modelo (solo la primera vez descarga ~1 GB)
device = "cuda" if torch.cuda.is_available() else "cpu"
model, preprocess = clip.load("ViT-B/32", device=device)

# Lista de posibles descripciones (puedes poner 1000 si quieres)
candidatos = [
    "un gato negro con gafas de sol",
    "un diagrama de Feynman con gluones",
    "un gráfico de dispersión con eje logarítmico",
    "una foto de un laboratorio de física de partículas",
    "una playa al atardecer",
    "una radiografía de pulmón",
    "un espectro de masas",
    "una pizarra llena de ecuaciones",
    "un perro corriendo en el parque",
    "una imagen de microscopía electrónica",
]

# Cargar imagen (desde URL o local)
url = "https://upload.wikimedia.org/wikipedia/commons/8/89/Feynman_Diagram_gluons.svg"  # cambia por tu imagen
# o desde archivo: Image.open("mi_imagen.jpg")
response = requests.get(url)
img = Image.open(BytesIO(response.content))

# Preprocesar
image_input = preprocess(img).unsqueeze(0).to(device)
text_tokens = clip.tokenize(candidatos).to(device)

with torch.no_grad():
    logits_per_image, _ = model(image_input, text_tokens)
    probs = logits_per_image.softmax(dim=-1).cpu().numpy()[0]

# Resultado
for texto, prob in zip(candidatos, probs):
    print(f"{prob*100:5.1f}% → {texto}")

# El ganador
mejor = candidatos[probs.argmax()]
print(f"\nDescripción más probable: {mejor}")
