import tkinter as tk
from tkinter import scrolledtext, messagebox
import threading
import ollama
import json
import os
from typing import List, Dict, Generator
import tempfile
import subprocess
import wave
import numpy as np
import sounddevice as sd

# --- CONFIGURACIÓN ---
MODEL_NAME = "gemma3"
PIPER_MODEL = "es_ES-carlfm-x-low.onnx"  # Cambia por el modelo Piper que tengas (ej. una voz española). Descárgalo de https://rhasspy.github.io/piper-samples/ o Hugging Face.
# Ejemplos: es_ES-m-ailabs-low.onnx, es_ES-sharerd-medium.onnx, etc.
# Asegúrate de tener el .onnx y su .json en la misma carpeta o especifica la ruta completa.

HABI_SYSTEM_PROMPT = """
# ... (el mismo prompt largo que tenías, lo omito aquí por brevedad pero déjalo igual)
"""

# ... (las clases SessionManager y OllamaClient permanecen iguales)

class HabiApp:
    def __init__(self, root: tk.Tk, ai_client: OllamaClient, storage: SessionManager):
        self.root = root
        self.ai = ai_client
        self.storage = storage
        self.history = []
        self.voice = None  # Se cargará la voz Piper al iniciar
        self.setup_ui()
        self.recover_session()
        self.load_piper_voice()

    def load_piper_voice(self):
        """Carga la voz de Piper una sola vez al iniciar la app."""
        try:
            from piper.voice import PiperVoice
            model_path = PIPER_MODEL  # Puedes poner ruta absoluta si no está en el directorio actual
            self.voice = PiperVoice.load(model_path)
            print(f"Voz Piper cargada: {PIPER_MODEL}")
        except Exception as e:
            messagebox.showerror("Error Piper TTS", f"No se pudo cargar la voz Piper:\n{e}\nInstala con: pip install piper-tts sounddevice numpy")
            self.voice = None

    def speak_text(self, text: str):
        """Reproduce el texto usando Piper TTS (síntesis + reproducción inmediata)."""
        if not self.voice or not text.strip():
            return

        threading.Thread(target=self._speak_worker, args=(text,), daemon=True).start()

    def _speak_worker(self, text: str):
        try:
            # Crear archivo WAV temporal
            with tempfile.NamedTemporaryFile(suffix=".wav", delete=False) as tmp_wav:
                wav_path = tmp_wav.name

            with wave.open(wav_path, "wb") as wav_file:
                # Configuración típica para la mayoría de modelos Piper (22050 Hz, mono, 16-bit)
                wav_file.setnchannels(1)
                wav_file.setsampwidth(2)  # 16-bit
                wav_file.setframerate(22050)
                self.voice.synthesize(text, wav_file)

            # Leer y reproducir con sounddevice
            data, fs = sd.read(wav_path, dtype='int16')
            sd.play(data, fs)
            sd.wait()  # Espera a que termine la reproducción

            # Limpiar archivo temporal
            os.unlink(wav_path)
        except Exception as e:
            print(f"Error en TTS: {e}")

    # ... (setup_ui permanece igual)

    def run_ai(self):
        full_reply = ""
        try:
            for chunk in self.ai.stream_chat(self.history):
                full_reply += chunk
                self.root.after(0, lambda c=chunk: self.append_text(c, None))
            self.history.append({"role": "assistant", "content": full_reply})
            self.root.after(0, lambda: self.append_text("\n\n", None))

            # --- NUEVO: Reproducir la respuesta completa en audio ---
            self.root.after(0, lambda: self.speak_text(full_reply))

            self.root.after(0, lambda: self.update_status(False))
        except Exception as e:
            error_msg = str(e)
            self.root.after(0, lambda m=error_msg: self.handle_error(m))
        finally:
            self.root.after(0, lambda: self.send_btn.config(state=tk.NORMAL))

    # ... (el resto de métodos permanecen iguales)

if __name__ == "__main__":
    # Dependencias adicionales necesarias:
    # pip install piper-tts sounddevice numpy
    app_root = tk.Tk()
    HabiApp(app_root, OllamaClient(model=MODEL_NAME), SessionManager())
    app_root.mainloop()
