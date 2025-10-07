#!/usr/bin/env python
import streamlit as st
import requests
import json
import os

# Configuración de Ollama (servidor local)
OLLAMA_URL = "http://localhost:11434/api/generate"
MODEL = "llama3.1"  # Cambia por tu modelo preferido

# Archivos externos para persistencia
HISTORIAL_JSON = "historial.json"
HISTORIAL_MD = "historial.md"

# Función para llamar a Ollama
def query_ollama(prompt, history=""):
    payload = {
        "model": MODEL,
        "prompt": f"{history}\nUsuario: {prompt}\nAsistente:",
        "stream": False
    }
    response = requests.post(OLLAMA_URL, json=payload)
    if response.status_code == 200:
        return response.json()["response"]
    else:
        return "Error: No se pudo conectar a Ollama."

# Cargar historial desde archivo JSON si existe
def cargar_historial():
    if os.path.exists(HISTORIAL_JSON):
        with open(HISTORIAL_JSON, "r") as f:
            return json.load(f)
    return []

# Guardar historial en JSON
def guardar_historial(messages):
    with open(HISTORIAL_JSON, "w") as f:
        json.dump(messages, f, ensure_ascii=False, indent=4)

# Exportar a Markdown (legible por humanos)
def exportar_a_markdown(messages):
    with open(HISTORIAL_MD, "w") as f:
        f.write("# Historial de Chat Cognitivo\n\n")
        for msg in messages:
            role = "Usuario" if msg["role"] == "user" else "Asistente"
            f.write(f"**{role}:** {msg['content']}\n\n---\n\n")

# Interfaz Streamlit
st.title("Tesela Cognitiva con Ollama")
st.write("Usa esta IA local como extensión de tu cognición. El historial se guarda en 'historial.json' (accesible a otras teselas) y se puede exportar a Markdown.")

# Cargar historial persistente
if "messages" not in st.session_state:
    st.session_state.messages = cargar_historial()

# Mostrar historial
for message in st.session_state.messages:
    with st.chat_message(message["role"]):
        st.markdown(message["content"])

# Entrada del usuario
if prompt := st.chat_input("¿Qué tarea cognitiva quieres distribuir?"):
    # Añadir mensaje del usuario
    st.session_state.messages.append({"role": "user", "content": prompt})
    with st.chat_message("user"):
        st.markdown(prompt)

    # Llamar a Ollama
    history = "\n".join([f"{msg['role']}: {msg['content']}" for msg in st.session_state.messages[:-1]])
    with st.chat_message("assistant"):
        response = query_ollama(prompt, history)
        st.markdown(response)
        st.session_state.messages.append({"role": "assistant", "content": response})

    # Guardar automáticamente en JSON (persistente y compartible)
    guardar_historial(st.session_state.messages)

# Botones para gestión
col1, col2 = st.columns(2)
with col1:
    if st.button("Resetear Tesela"):
        st.session_state.messages = []
        guardar_historial([])  # Limpiar archivo
        if os.path.exists(HISTORIAL_MD):
            os.remove(HISTORIAL_MD)
        st.rerun()
with col2:
    if st.button("Exportar a Markdown"):
        exportar_a_markdown(st.session_state.messages)
        st.success(f"Historial exportado a {HISTORIAL_MD} (legible por ti).")
