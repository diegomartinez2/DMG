import tkinter as tk
from tkinter import scrolledtext, messagebox
import threading
import ollama
import json
import os
from typing import List, Dict, Generator

# --- CONFIGURACIÓN ---
MODEL_NAME = "gemma3"

# Prompt del sistema especializado en programación
HABI_SYSTEM_PROMPT = """Eres ‘Habi’, una experta asistente de programación.
Tu objetivo es ayudar a modificar, depurar y optimizar el código que el usuario te proporciona en el editor lateral.
* No utilices emojis.
* Responde en el mismo idioma que el usuario.
* Sé directa: si el usuario pide un cambio, muestra el código corregido o explica el error.
* Utiliza siempre el "CONTEXTO DEL CÓDIGO ACTUAL" como referencia principal."""

class SessionManager:
    def __init__(self, filename: str = "habi_code_history.json"):
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

    def stream_chat(self, history: List[Dict[str, str]], code_content: str) -> Generator[str, None, None]:
        # Evitamos el error de "unterminated f-string" pasando el código
        # como un bloque de texto independiente antes de las f-strings.
        system_instruction = f"{HABI_SYSTEM_PROMPT}\n\nCONTEXTO DEL CÓDIGO ACTUAL:\n{code_content}"

        # Construimos la lista de mensajes
        messages = [{'role': 'system', 'content': system_instruction}] + history

        try:
            response = ollama.chat(model=self.model, messages=messages, stream=True)
            for chunk in response:
                yield chunk['message']['content']
        except Exception as e:
            raise ConnectionError(f"Ollama Error: {str(e)}")

class HabiApp:
    def __init__(self, root: tk.Tk):
        self.root = root
        self.ai = OllamaClient()
        self.storage = SessionManager()
        self.history = []

        self.setup_ui()
        self.load_history()

    def setup_ui(self):
        self.root.title(f"Habi - Programación (Gemma 3)")
        self.root.geometry("1100x700")
        self.root.configure(bg="#1e1e1e")

        # PanedWindow para dividir Chat y Código
        self.paned = tk.PanedWindow(self.root, orient=tk.HORIZONTAL, bg="#333", sashwidth=4)
        self.paned.pack(fill=tk.BOTH, expand=True)

        # LADO IZQUIERDO: Chat
        self.chat_frame = tk.Frame(self.paned, bg="#1e1e1e")
        self.paned.add(self.chat_frame, width=450)

        self.display = scrolledtext.ScrolledText(self.chat_frame, bg="#1e1e1e", fg="#d4d4d4",
                                                font=("Consolas", 10), state=tk.DISABLED)
        self.display.pack(padx=10, pady=10, fill=tk.BOTH, expand=True)
        self.display.tag_config("user", foreground="#569cd6", font=("Consolas", 10, "bold"))

        self.input_text = tk.Text(self.chat_frame, height=3, bg="#2d2d2d", fg="white", font=("Consolas", 10))
        self.input_text.pack(padx=10, pady=(0, 10), fill=tk.X)
        self.input_text.bind("<Return>", self.handle_enter)

        self.send_btn = tk.Button(self.chat_frame, text="Preguntar a Habi", command=self.on_send,
                                 bg="#0e639c", fg="white", relief=tk.FLAT)
        self.send_btn.pack(padx=10, pady=(0, 10), fill=tk.X)

        # LADO DERECHO: Editor de Código
        self.code_frame = tk.Frame(self.paned, bg="#1e1e1e")
        self.paned.add(self.code_frame, width=650)

        tk.Label(self.code_frame, text="SCRIPT A MODIFICAR", bg="#1e1e1e", fg="#808080", font=("Arial", 9)).pack()

        self.code_editor = scrolledtext.ScrolledText(self.code_frame, bg="#1e1e1e", fg="#ce9178",
                                                    font=("Consolas", 11), insertbackground="white")
        self.code_editor.pack(padx=10, pady=10, fill=tk.BOTH, expand=True)
        self.code_editor.insert("1.0", "# Pega aquí el código que quieres que Habi analice...\n")

    def handle_enter(self, event):
        if not (event.state & 0x1): # Si no es Shift
            self.on_send()
            return "break"

    def append_text(self, content, tag=None):
        self.display.config(state=tk.NORMAL)
        self.display.insert(tk.END, content, tag)
        self.display.see(tk.END)
        self.display.config(state=tk.DISABLED)

    def on_send(self):
        prompt = self.input_text.get("1.0", tk.END).strip()
        if not prompt: return

        # Obtenemos el código actual para que la IA lo lea en este turno
        current_code = self.code_editor.get("1.0", tk.END).strip()

        self.input_text.delete("1.0", tk.END)
        self.history.append({"role": "user", "content": prompt})
        self.append_text(f"Tú: {prompt}\n\nHabi: ", "user")

        self.send_btn.config(state=tk.DISABLED)
        threading.Thread(target=self.run_ai, args=(current_code,), daemon=True).start()

    def run_ai(self, code_context):
        full_reply = ""
        try:
            for chunk in self.ai.stream_chat(self.history, code_context):
                full_reply += chunk
                self.root.after(0, lambda c=chunk: self.append_text(c))

            self.history.append({"role": "assistant", "content": full_reply})
            self.root.after(0, lambda: self.append_text("\n\n"))
            self.storage.save(self.history)
        except Exception as e:
            self.root.after(0, lambda: messagebox.showerror("Error", str(e)))
        finally:
            self.root.after(0, lambda: self.send_btn.config(state=tk.NORMAL))

    def load_history(self):
        self.history = self.storage.load()
        for msg in self.history:
            tag = "user" if msg['role'] == 'user' else None
            label = "Tú: " if msg['role'] == 'user' else "Habi: "
            self.append_text(f"{label}{msg['content']}\n\n", tag)

if __name__ == "__main__":
    app_root = tk.Tk()
    HabiApp(app_root)
    app_root.mainloop()
