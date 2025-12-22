import tkinter as tk
from tkinter import scrolledtext, messagebox
import threading
import ollama
import json
import os
import wave
import subprocess
from typing import List, Dict, Generator

# Intentar importar Piper
try:
    from piper.voice import PiperVoice
except ImportError:
    print("Error: La librería 'piper-tts' no está instalada.")

# --- CONFIGURACIÓN ---
MODEL_NAME = "gemma3"
PIPER_MODEL_PATH = os.path.expanduser("~/piper_voices/es_ES-sharvard-medium.onnx")

HABI_SYSTEM_PROMPT = """Eres ‘Habi’, una 'musa AI' brillante y amigable.
Tu propósito es ayudar al usuario a comprender y aprender.
* No utilices emojis.
* Responde en el mismo idioma que el usuario.
* Nunca narres lo que vas a hacer, simplemente hazlo."""

class PiperEngine:
    """Clase para manejar la síntesis de voz local con Piper."""
    def __init__(self, model_path: str):
        self.model_path = model_path
        self.voice = None
        self._load_voice()

    def _load_voice(self):
        if os.path.exists(self.model_path):
            try:
                self.voice = PiperVoice.load(self.model_path)
                print(f"Voz de Piper cargada correctamente: {self.model_path}")
            except Exception as e:
                print(f"Error cargando voz de Piper: {e}")
        else:
            print(f"Aviso: No se encontró el modelo en {self.model_path}")

    def speak(self, text: str):
        """Sintetiza texto asegurando la compatibilidad con objetos AudioChunk."""
        if not self.voice:
            return

        def task():
            output_file = "temp_voice.wav"
            try:
                # Piper puede escribir directamente en un archivo wave abierto
                with wave.open(output_file, "wb") as wav_file:
                    # Configuración estandar para Piper
                    wav_file.setnchannels(1)
                    wav_file.setsampwidth(2) # 16-bit
                    wav_file.setframerate(self.voice.config.sample_rate)

                    # Iteramos sobre la síntesis
                    for chunk in self.voice.synthesize(text):
                        # Intentamos obtener los bytes de varias formas posibles según la versión
                        if isinstance(chunk, bytes):
                            wav_file.writeframes(chunk)
                        elif hasattr(chunk, 'audio_data'): # Caso común en versiones recientes
                            wav_file.writeframes(chunk.audio_data)
                        elif hasattr(chunk, 'audio'): # Otra posible variante
                            wav_file.writeframes(chunk.audio)
                        else:
                            # Si es un objeto AudioChunk pero no sabemos cómo leerlo,
                            # forzamos la conversión a bytes si Piper lo permite
                            try:
                                wav_file.writeframes(bytes(chunk))
                            except:
                                print(f"No se pudo procesar el fragmento de audio: {type(chunk)}")

                # Reproducir usando aplay (estándar en Linux)
                if os.path.exists(output_file) and os.path.getsize(output_file) > 44:
                    subprocess.run(["aplay", output_file], check=True, stderr=subprocess.DEVNULL)

            except Exception as e:
                print(f"Error en la síntesis de voz: {e}")
            finally:
                if os.path.exists(output_file):
                    try: os.remove(output_file)
                    except: pass

        threading.Thread(target=task, daemon=True).start()

class SessionManager:
    def __init__(self, filename: str = "habi_history.json"):
        self.filename = filename

    def save(self, history: List[Dict[str, str]]):
        try:
            with open(self.filename, 'w', encoding='utf-8') as f:
                json.dump(history, f, ensure_ascii=False, indent=4)
        except Exception as e:
            print(f"Error al guardar historial: {e}")

    def load(self) -> List[Dict[str, str]]:
        if os.path.exists(self.filename):
            try:
                with open(self.filename, 'r', encoding='utf-8') as f:
                    return json.load(f)
            except:
                return []
        return []

class OllamaClient:
    def __init__(self, model: str = MODEL_NAME):
        self.model = model

    def stream_chat(self, history: List[Dict[str, str]]) -> Generator[str, None, None]:
        messages = [{'role': 'system', 'content': HABI_SYSTEM_PROMPT}] + history
        try:
            response = ollama.chat(model=self.model, messages=messages, stream=True)
            for chunk in response:
                yield chunk['message']['content']
        except Exception as e:
            raise ConnectionError(f"Ollama no responde: {e}")

class HabiApp:
    def __init__(self, root: tk.Tk, ai_client: OllamaClient, storage: SessionManager, piper: PiperEngine):
        self.root = root
        self.ai = ai_client
        self.storage = storage
        self.piper = piper
        self.history = []

        self.setup_ui()
        self.recover_session()

    def setup_ui(self):
        self.root.title(f"Habi AI - {MODEL_NAME}")
        self.root.geometry("800x800")
        self.root.configure(bg="#0f172a")

        self.display = scrolledtext.ScrolledText(self.root, bg="#1e293b", fg="#f1f5f9", font=("Segoe UI", 11), state=tk.DISABLED)
        self.display.pack(padx=20, pady=(20, 10), fill=tk.BOTH, expand=True)
        self.display.tag_config("user", foreground="#38bdf8", font=("Segoe UI", 10, "bold"))
        self.display.tag_config("error", foreground="#f87171")

        self.input_text = tk.Text(self.root, height=3, bg="#1e293b", fg="white", font=("Segoe UI", 11), insertbackground="white")
        self.input_text.pack(padx=20, pady=(0, 10), fill=tk.X)
        self.input_text.bind("<Return>", self.handle_return)

        self.status_frame = tk.Frame(self.root, bg="#0f172a")
        self.status_frame.pack(padx=20, pady=(0, 10), fill=tk.X)

        self.status_light = tk.Canvas(self.status_frame, width=15, height=15, bg="#0f172a", highlightthickness=0)
        self.status_light.pack(side=tk.LEFT)
        self.light_id = self.status_light.create_oval(2, 2, 13, 13, fill="#ef4444")

        self.status_label = tk.Label(self.status_frame, text="Esperando...", bg="#0f172a", fg="#94a3b8", font=("Segoe UI", 9))
        self.status_label.pack(side=tk.LEFT, padx=5)

        btn_frame = tk.Frame(self.root, bg="#0f172a")
        btn_frame.pack(padx=20, pady=(0, 20), fill=tk.X)

        self.send_btn = tk.Button(btn_frame, text="Enviar", command=self.on_send, bg="#0284c7", fg="white", relief=tk.FLAT, font=("Segoe UI", 10, "bold"))
        self.send_btn.pack(side=tk.LEFT, expand=True, fill=tk.X)

        tk.Button(btn_frame, text="Guardar", command=self.on_save, bg="#059669", fg="white", relief=tk.FLAT).pack(side=tk.RIGHT, expand=True, fill=tk.X, padx=(10, 0))

    def update_status(self, active: bool, message: str = ""):
        color = "#22c55e" if active else "#ef4444"
        self.status_light.itemconfig(self.light_id, fill=color)
        self.status_label.config(text=message if message else ("Pensando..." if active else "Listo"))

    def handle_return(self, event):
        if not (event.state & 0x1):
            self.on_send()
            return "break"

    def recover_session(self):
        data = self.storage.load()
        if data:
            self.history = data
            for msg in self.history:
                tag = "user" if msg['role'] == 'user' else None
                self.append_text(f"{'Tú' if msg['role'] == 'user' else 'Habi'}: {msg['content']}\n\n", tag)

    def append_text(self, content, tag):
        self.display.config(state=tk.NORMAL)
        self.display.insert(tk.END, content, tag)
        self.display.see(tk.END)
        self.display.config(state=tk.DISABLED)

    def on_save(self):
        self.storage.save(self.history)
        messagebox.showinfo("Sistema", "Historial guardado.")

    def on_send(self):
        text = self.input_text.get("1.0", tk.END).strip()
        if not text: return

        self.input_text.delete("1.0", tk.END)
        self.history.append({"role": "user", "content": text})
        self.append_text(f"Tú: {text}\n\nHabi: ", "user")

        self.send_btn.config(state=tk.DISABLED)
        self.update_status(True)
        threading.Thread(target=self.run_ai, daemon=True).start()

    def run_ai(self):
        full_reply = ""
        try:
            for chunk in self.ai.stream_chat(self.history):
                full_reply += chunk
                self.root.after(0, lambda c=chunk: self.append_text(c, None))

            self.history.append({"role": "assistant", "content": full_reply})
            self.root.after(0, lambda: self.append_text("\n\n", None))
            self.piper.speak(full_reply)
            self.root.after(0, lambda: self.update_status(False))
        except Exception as e:
            self.root.after(0, lambda: self.handle_error(str(e)))
        finally:
            self.root.after(0, lambda: self.send_btn.config(state=tk.NORMAL))

    def handle_error(self, message):
        self.update_status(False, "Error")
        self.append_text(f"\n[ERROR]: {message}\n\n", "error")

if __name__ == "__main__":
    root = tk.Tk()
    engine = PiperEngine(PIPER_MODEL_PATH)
    client = OllamaClient()
    storage = SessionManager()
    app = HabiApp(root, client, storage, engine)
    root.mainloop()
