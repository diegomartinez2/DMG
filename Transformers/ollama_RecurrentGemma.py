import tkinter as tk
from tkinter import scrolledtext, messagebox
import threading
import ollama
import json
import os
import wave
import subprocess
from typing import List, Dict, Generator
from piper.voice import PiperVoice

# --- CONFIGURACIÓN ---
MODEL_NAME = "gemma3"
# Ruta al modelo de Piper (Ajusta según tu instalación)
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
            except Exception as e:
                print(f"Error cargando voz de Piper: {e}")

    def speak(self, text: str):
        """Sintetiza texto y lo reproduce usando aplay."""
        if not self.voice:
            print("Voz no cargada. Piper no está disponible.")
            return

        def task():
            output_file = "temp_voice.wav"
            try:
                with wave.open(output_file, "wb") as wav_file:
                    self.voice.synthesize(text, wav_file)
                # Reproducir usando aplay (estándar en Ubuntu)
                subprocess.run(["aplay", output_file], stderr=subprocess.DEVNULL)
            except Exception as e:
                print(f"Error en reproducción de voz: {e}")
            finally:
                if os.path.exists(output_file):
                    try: os.remove(output_file)
                    except: pass

        # Ejecutar en hilo separado para no bloquear la UI o el stream
        threading.Thread(target=task, daemon=True).start()

class SessionManager:
    def __init__(self, filename: str = "habi_history.json"):
        self.filename = filename

    def save(self, history: List[Dict[str, str]]):
        try:
            with open(self.filename, 'w', encoding='utf-8') as f:
                json.dump(history, f, ensure_ascii=False, indent=4)
        except Exception as e:
            print(f"Error al guardar: {e}")

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
            raise ConnectionError("No se pudo conectar con Ollama. ¿Está el servidor encendido?") from e

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
        self.root.title(f"Habi AI con Voz - {MODEL_NAME}")
        self.root.geometry("800x750")
        self.root.configure(bg="#0f172a")

        self.display = scrolledtext.ScrolledText(self.root, bg="#1e293b", fg="#f1f5f9", font=("Segoe UI", 11), state=tk.DISABLED)
        self.display.pack(padx=20, pady=(20, 10), fill=tk.BOTH, expand=True)
        self.display.tag_config("user", foreground="#38bdf8", font=("Segoe UI", 10, "bold"))
        self.display.tag_config("error", foreground="#f87171", font=("Segoe UI", 10, "italic"))

        self.input_text = tk.Text(self.root, height=3, bg="#1e293b", fg="white", font=("Segoe UI", 11), insertbackground="white")
        self.input_text.pack(padx=20, pady=(0, 10), fill=tk.X)
        self.input_text.bind("<Return>", self.handle_return)

        self.status_frame = tk.Frame(self.root, bg="#0f172a")
        self.status_frame.pack(padx=20, pady=(0, 10), fill=tk.X)

        self.status_light = tk.Canvas(self.status_frame, width=20, height=20, bg="#0f172a", highlightthickness=0)
        self.status_light.pack(side=tk.LEFT)
        self.light_id = self.status_light.create_oval(4, 4, 16, 16, fill="#ef4444")

        self.status_label = tk.Label(self.status_frame, text="Listo / En reposo", bg="#0f172a", fg="#94a3b8", font=("Segoe UI", 9))
        self.status_label.pack(side=tk.LEFT, padx=5)

        btn_frame = tk.Frame(self.root, bg="#0f172a")
        btn_frame.pack(padx=20, pady=(0, 20), fill=tk.X)

        self.send_btn = tk.Button(btn_frame, text="Preguntar", command=self.on_send, bg="#0284c7", fg="white", font=("Segoe UI", 10, "bold"), relief=tk.FLAT)
        self.send_btn.pack(side=tk.LEFT, expand=True, fill=tk.X)

        tk.Button(btn_frame, text="Guardar Sesión", command=self.on_save, bg="#059669", fg="white", font=("Segoe UI", 10, "bold"), relief=tk.FLAT).pack(side=tk.RIGHT, expand=True, fill=tk.X, padx=(10, 0))

    def update_status(self, active: bool, message: str = ""):
        color = "#22c55e" if active else "#ef4444"
        default_msg = "Habi está pensando y hablando..." if active else "Listo / En reposo"
        self.status_light.itemconfig(self.light_id, fill=color)
        self.status_label.config(text=message if message else default_msg)

    def handle_return(self, event):
        if not (event.state & 0x1):
            self.on_send()
            return "break"

    def recover_session(self):
        data = self.storage.load()
        if data:
            self.history = data
            self.append_text("--- Sesión previa cargada ---\n\n", None)
            for msg in self.history:
                tag = "user" if msg['role'] == 'user' else None
                label = "Tú: " if msg['role'] == 'user' else "Habi: "
                self.append_text(f"{label}{msg['content']}\n\n", tag)

    def append_text(self, content, tag):
        self.display.config(state=tk.NORMAL)
        self.display.insert(tk.END, content, tag)
        self.display.see(tk.END)
        self.display.config(state=tk.DISABLED)

    def on_save(self):
        self.storage.save(self.history)
        messagebox.showinfo("Habi", "Historial guardado correctamente.")

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

            # Una vez terminada la respuesta de texto, activar la voz
            self.history.append({"role": "assistant", "content": full_reply})
            self.root.after(0, lambda: self.append_text("\n\n", None))

            # Hablar la respuesta completa
            self.piper.speak(full_reply)

            self.root.after(0, lambda: self.update_status(False))
        except ConnectionError as e:
            self.root.after(0, lambda: self.handle_error(str(e)))
        except Exception as e:
            self.root.after(0, lambda: self.handle_error(f"Error inesperado: {str(e)}"))
        finally:
            self.root.after(0, lambda: self.send_btn.config(state=tk.NORMAL))

    def handle_error(self, message):
        self.update_status(False, "Error")
        self.append_text(f"\n[SISTEMA]: {message}\n\n", "error")

if __name__ == "__main__":
    app_root = tk.Tk()
    piper_engine = PiperEngine(PIPER_MODEL_PATH)
    app_root.mainloop()
    HabiApp(app_root, OllamaClient(model=MODEL_NAME), SessionManager(), piper_engine)
    app_root.mainloop()
