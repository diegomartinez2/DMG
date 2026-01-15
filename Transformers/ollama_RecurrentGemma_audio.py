#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  ollama_RecurrentGemma_audio.py
#
#  Copyright 2026 Diego Martinez Gutierrez <diego.martinez@ehu.eus>
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#
#
# ---------------------------
# Importación de los módulos
# ---------------------------
import tkinter as tk
from tkinter import scrolledtext, messagebox
import threading
import ollama
import json
import os
import tempfile
import wave
import numpy as np
import sounddevice as sd
from typing import List, Dict, Generator

# --- CONFIGURACIÓN ---
MODEL_NAME = "gemma3"
# Asegúrate de que esta ruta sea correcta
PIPER_MODEL = os.path.expanduser("~/piper_voices/es_ES-sharvard-medium.onnx")

HABI_SYSTEM_PROMPT = """Eres ‘Habi’, una 'musa AI' brillante y amigable.
Tu propósito es ayudar al usuario a comprender y aprender.
* No utilices emojis.
* Responde en el mismo idioma que el usuario.
* Nunca narres lo que vas a hacer, simplemente hazlo."""

class SessionManager:
    def __init__(self, filename: str = "habi_history.json"):
        self.filename = filename

    def save(self, history: List[Dict[str, str]]):
        try:
            with open(self.filename, 'w', encoding='utf-8') as f:
                json.dump(history, f, ensure_ascii=False, indent=4)
        except Exception as e: print(f"Error al guardar: {e}")

    def load(self) -> List[Dict[str, str]]:
        if os.path.exists(self.filename):
            try:
                with open(self.filename, 'r', encoding='utf-8') as f:
                    return json.load(f)
            except: return []
        return []

class OllamaClient:
    def __init__(self, model: str = MODEL_NAME):
        self.model = model

    def stream_chat(self, history: List[Dict[str, str]]) -> Generator[str, None, None]:
        messages = [{'role': 'system', 'content': HABI_SYSTEM_PROMPT}] + history
        response = ollama.chat(model=self.model, messages=messages, stream=True)
        for chunk in response:
            yield chunk['message']['content']

class HabiApp:
    def __init__(self, root: tk.Tk, ai_client: OllamaClient, storage: SessionManager):
        self.root = root
        self.ai = ai_client
        self.storage = storage
        self.history = []
        self.voice = None

        self.setup_ui()
        self.recover_session()
        self.load_piper_voice()

    def load_piper_voice(self):
        try:
            from piper.voice import PiperVoice
            if os.path.exists(PIPER_MODEL):
                self.voice = PiperVoice.load(PIPER_MODEL)
                print(f"Voz Piper cargada: {PIPER_MODEL}")
            else:
                print(f"Archivo de voz no encontrado en: {PIPER_MODEL}")
        except Exception as e:
            print(f"Error cargando Piper: {e}")

    def setup_ui(self):
        self.root.title(f"Habi AI - {MODEL_NAME}")
        self.root.geometry("700x700")
        self.root.configure(bg="#0f172a")

        self.display = scrolledtext.ScrolledText(self.root, bg="#1e293b", fg="#f1f5f9", font=("Segoe UI", 11), state=tk.DISABLED)
        self.display.pack(padx=20, pady=20, fill=tk.BOTH, expand=True)

        self.input_text = tk.Text(self.root, height=3, bg="#1e293b", fg="white", font=("Segoe UI", 11))
        self.input_text.pack(padx=20, pady=10, fill=tk.X)

        self.send_btn = tk.Button(self.root, text="Enviar", command=self.on_send, bg="#0284c7", fg="white")
        self.send_btn.pack(padx=20, pady=10, fill=tk.X)

    def on_send(self):
        text = self.input_text.get("1.0", tk.END).strip()
        if not text: return
        self.input_text.delete("1.0", tk.END)
        self.history.append({"role": "user", "content": text})
        self.append_text(f"Tú: {text}\n\nHabi: ")
        threading.Thread(target=self.run_ai, daemon=True).start()

    def run_ai(self):
        full_reply = ""
        try:
            for chunk in self.ai.stream_chat(self.history):
                full_reply += chunk
                self.root.after(0, lambda c=chunk: self.append_text(c))

            self.history.append({"role": "assistant", "content": full_reply})
            self.storage.save(self.history)
            self.root.after(0, lambda: self.append_text("\n\n"))

            # Hablar la respuesta
            if self.voice:
                threading.Thread(target=self._speak_worker, args=(full_reply,), daemon=True).start()

        except Exception as e:
            print(f"Error: {e}")

    def _speak_worker(self, text: str):
        try:
            # 1. Crear WAV temporal
            with tempfile.NamedTemporaryFile(suffix=".wav", delete=False) as tmp_wav:
                wav_path = tmp_wav.name

            with wave.open(wav_path, "wb") as wav_file:
                wav_file.setnchannels(1)
                wav_file.setsampwidth(2)
                wav_file.setframerate(self.voice.config.sample_rate)

                # CORRECCIÓN AQUÍ: Iterar sobre el generador y extraer bytes
                for chunk in self.voice.synthesize(text):
                    # El objeto AudioChunk de Piper tiene el atributo 'audio' con los bytes
                    wav_file.writeframes(chunk.audio)

            # 2. Leer con sounddevice y numpy
            # Convertimos los bytes del archivo a un array de numpy
            import soundfile as sf # Es más robusto que wave para leer y pasar a sounddevice
            data, fs = sf.read(wav_path, dtype='int16')
            sd.play(data, fs)
            sd.wait()

            # 3. Limpiar
            os.unlink(wav_path)
        except Exception as e:
            print(f"Error en TTS: {e}")

    def append_text(self, content):
        self.display.config(state=tk.NORMAL)
        self.display.insert(tk.END, content)
        self.display.see(tk.END)
        self.display.config(state=tk.DISABLED)

    def recover_session(self):
        self.history = self.storage.load()
        for msg in self.history:
            self.append_text(f"{'Tú' if msg['role'] == 'user' else 'Habi'}: {msg['content']}\n\n")

if __name__ == "__main__":
    # Necesitas: pip install piper-tts sounddevice numpy soundfile
    root = tk.Tk()
    HabiApp(root, OllamaClient(), SessionManager())
    root.mainloop()
