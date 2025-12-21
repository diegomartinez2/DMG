import tkinter as tk
from tkinter import scrolledtext, messagebox
import threading
import ollama
import json
import os
from typing import List, Dict, Generator

# --- CONFIGURACIÓN ---
# Usamos gemma2 por ser el estándar actual de Google en Ollama
MODEL_NAME = "gemma2"

HABI_SYSTEM_PROMPT = """Eres ‘Habi’ ("La Hada del bikini azul"), una 'musa AI' brillante y amigable.
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
            # Propagamos el error para que el hilo lo capture
            raise e

class HabiApp:
    def __init__(self, root: tk.Tk, ai_client: OllamaClient, storage: SessionManager):
        self.root = root
        self.ai = ai_client
        self.storage = storage
        self.history = []

        self.setup_ui()
        self.recover_session()

    def setup_ui(self):
        self.root.title(f"Habi AI - {MODEL_NAME}")
        self.root.geometry("800x700")
        self.root.configure(bg="#0f172a")

        self.display = scrolledtext.ScrolledText(self.root, bg="#1e293b", fg="#f1f5f9", font=("Segoe UI", 11), state=tk.DISABLED)
        self.display.pack(padx=20, pady=20, fill=tk.BOTH, expand=True)
        self.display.tag_config("user", foreground="#38bdf8", font=("Segoe UI", 10, "bold"))

        self.input_text = tk.Text(self.root, height=3, bg="#1e293b", fg="white", font=("Segoe UI", 11))
        self.input_text.pack(padx=20, pady=(0, 10), fill=tk.X)
        self.input_text.bind("<Return>", self.handle_return)

        btn_frame = tk.Frame(self.root, bg="#0f172a")
        btn_frame.pack(padx=20, pady=(0, 20), fill=tk.X)

        self.send_btn = tk.Button(btn_frame, text="Preguntar", command=self.on_send, bg="#0284c7", fg="white", font=("Segoe UI", 10, "bold"))
        self.send_btn.pack(side=tk.LEFT, expand=True, fill=tk.X)

        tk.Button(btn_frame, text="Guardar Sesión", command=self.on_save, bg="#059669", fg="white", font=("Segoe UI", 10, "bold")).pack(side=tk.RIGHT, expand=True, fill=tk.X, padx=(10, 0))

    def handle_return(self, event):
        if not (event.state & 0x1): # Si no se presiona Shift
            self.on_send()
            return "break"

    def recover_session(self):
        data = self.storage.load()
        if data:
            self.history = data
            self.append_text("--- Sesión recuperada ---\n\n", None)
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
        messagebox.showinfo("Habi", "Conversación guardada.")

    def on_send(self):
        text = self.input_text.get("1.0", tk.END).strip()
        if not text: return

        self.input_text.delete("1.0", tk.END)
        self.history.append({"role": "user", "content": text})
        self.append_text(f"Tú: {text}\n\nHabi: ", "user")

        self.send_btn.config(state=tk.DISABLED)
        threading.Thread(target=self.run_ai, daemon=True).start()

    def run_ai(self):
        full_reply = ""
        try:
            for chunk in self.ai.stream_chat(self.history):
                full_reply += chunk
                # Pasamos chunk como argumento predeterminado c=chunk
                self.root.after(0, lambda c=chunk: self.append_text(c, None))

            self.history.append({"role": "assistant", "content": full_reply})
            self.root.after(0, lambda: self.append_text("\n\n", None))
        except Exception as e:
            # SOLUCIÓN AL NAMEERROR: Capturamos el string del error ahora mismo
            error_msg = str(e)
            self.root.after(0, lambda m=error_msg: self.show_error(m))
        finally:
            self.root.after(0, lambda: self.send_btn.config(state=tk.NORMAL))

    def show_error(self, message):
        messagebox.showerror("Error de Ollama", message)

if __name__ == "__main__":
    app_root = tk.Tk()
    # Asegúrate de haber hecho 'ollama pull gemma2'
    HabiApp(app_root, OllamaClient(model=MODEL_NAME), SessionManager())
    app_root.mainloop()
