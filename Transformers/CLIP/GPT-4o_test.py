# pip install openai pillow
from openai import OpenAI
from PIL import Image
import base64

#client = OpenAI(api_key="tu_api_key_aquí")  # o usa Groq gratis con clave de groq.com
client = OpenAI(
    api_key="gsk_tu_clave_aquí",  # ¡Pégala aquí!
    base_url="https://api.groq.com/openai/v1"  # Esto es clave para Groq
)
def encode_image(image_path):
    with open(image_path, "rb") as f:
        return base64.b64encode(f.read()).decode('utf-8')

img_path = "mi_grafico_fisica.jpg"
base64_image = encode_image(img_path)

response = client.chat.completions.create(
    model="gpt-4o",
    messages=[
        {"role": "user", "content": [
            {"type": "text", "text": "Extrae todos los datos numéricos de este gráfico y dame una tabla markdown. Luego explica qué muestra el experimento."},
            {"type": "image_url", "image_url": {"url": f"data:image/jpeg;base64,{base64_image}"}}
        ]}
    ],
    max_tokens=1000
)

print(response.choices[0].message.content)
