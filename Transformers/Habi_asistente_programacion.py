import tkinter as tk
from tkinter import scrolledtext, messagebox, font
import threading
import ollama
import json
import os
from typing import List, Dict, Generator

# --- CONFIGURACIÓN ---
MODEL_NAME = "gemma3"

HABI_SYSTEM_PROMPT = """Eres ‘Habi’, una musa AI y experta programadora.
Tu propósito es ayudar al usuario a escribir, depurar y optimizar código siguiendo los principios SOLID y KISS.
* No utilices emojis.
* Responde en el mismo idioma que el usuario.
* Nunca narres lo que vas a hacer, simplemente hazlo.
* Cuando el usuario te pase un código, analízalo y propón mejoras o realiza las modificaciones solicitadas directamente."""

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

    def stream_chat(self, history: List[Dict[str, str]], code_context: str) -> Generator[str, None, None]:
        # Construimos el mensaje incluyendo el contexto del código actual
        context_msg = f"CONTEXTO DEL CÓDIGO ACTUAL:\n

        # Insertamos el contexto del código justo después del system prompt para que siempre esté presente
        messages = [
            {'role': 'system', 'content': HABI_SYSTEM_PROMPT},
            {'role': 'system', 'content': context_msg}
        ] + history

        try:
            response = ollama.chat(model=self.model, messages=messages, stream=True)
            for chunk in response:
                yield chunk['message']['content']
        except Exception as e:
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
        self.root.title(f"Habi AI - Code Assistant ({MODEL_NAME})")
        self.root.geometry("1200x800")
        self.root.configure(bg="#0f172a")

        # Layout Principal: Ventana dividida (PanedWindow)
        self.paned = tk.PanedWindow(self.root, orient=tk.HORIZONTAL, bg="#1e293b", sashwidth=4)
        self.paned.pack(fill=tk.BOTH, expand=True)

        # --- LADO IZQUIERDO: CHAT ---
        self.chat_frame = tk.Frame(self.paned, bg="#0f172a")
        self.paned.add(self.chat_frame, width=500)

        self.display = scrolledtext.ScrolledText(self.chat_frame, bg="#0f172a", fg="#f1f5f9",
                                                font=("Consolas", 10), state=tk.DISABLED, borderwidth=0)
        self.display.pack(padx=10, pady=10, fill=tk.BOTH, expand=True)
        self.display.tag_config("user", foreground="#38bdf8", font=("Consolas", 10, "bold"))

        self.input_text = tk.Text(self.chat_frame, height=4, bg="#1e293b", fg="white",
                                 font=("Consolas", 10), insertbackground="white")
        self.input_text.pack(padx=10, pady=(0, 10), fill=tk.X)
        self.input_text.bind("<Return>", self.handle_return)

        # Botones Chat
        btn_chat_frame = tk.Frame(self.chat_frame, bg="#0f172a")
        btn_chat_frame.pack(padx=10, pady=(0, 10), fill=tk.X)

        self.send_btn = tk.Button(btn_chat_frame, text="Enviar Prompt", command=self.on_send,
                                 bg="#0284c7", fg="white", font=("Segoe UI", 9, "bold"), relief=tk.FLAT)
        self.send_btn.pack(side=tk.LEFT, expand=True, fill=tk.X)

        # --- LADO DERECHO: EDITOR DE CÓDIGO ---
        self.code_frame = tk.Frame(self.paned, bg="#1e293b")
        self.paned.add(self.code_frame, width=700)

        tk.Label(self.code_frame, text="EDITOR DE CÓDIGO (Contexto para Habi)", bg="#1e293b",
                 fg="#94a3b8", font=("Segoe UI", 9, "bold")).pack(pady=5)

        self.code_editor = scrolledtext.ScrolledText(self.code_frame, bg="#000000", fg="#adff2f",
                                                    font=("Consolas", 11), insertbackground="white",
                                                    undo=True, borderwidth=0)
        self.code_editor.pack(padx=10, pady=(0, 10), fill=tk.BOTH, expand=True)
        # Código de ejemplo inicial
        self.code_editor.insert("1.0", "# Pega aquí el código que quieres modificar...\ndef ejemplo():\n    print('Hola mundo')\n")

    def handle_return(self, event):
        if not (event.state & 0x1): # Si no es Shift+Enter
            self.on_send()
            return "break"

    def append_text(self, content, tag):
        self.display.config(state=tk.NORMAL)
        self.display.insert(tk.END, content, tag)
        self.display.see(tk.END)
        self.display.config(state=tk.DISABLED)

    def on_send(self):
        prompt = self.input_text.get("1.0", tk.END).strip()
        if not prompt: return

        # Obtener el código actual del editor para enviarlo como contexto
        current_code = self.code_editor.get("1.0", tk.END).strip()

        self.input_text.delete("1.0", tk.END)
        self.history.append({"role": "user", "content": prompt})
        self.append_text(f"Tú: {prompt}\n\nHabi: ", "user")

        self.send_btn.config(state=tk.DISABLED)

        # Pasamos el contexto del código a la función de ejecución
        threading.Thread(target=self.run_ai, args=(current_code,), daemon=True).start()

    def run_ai(self, code_context):
        full_reply = ""
        try:
            for chunk in self.ai.stream_chat(self.history, code_context):
                full_reply += chunk
                self.root.after(0, lambda c=chunk: self.append_text(c, None))

            self.history.append({"role": "assistant", "content": full_reply})
            self.root.after(0, lambda: self.append_text("\n\n---\n\n", None))
        except Exception as e:
            self.root.after(0, lambda m=str(e): messagebox.showerror("Error", m))
        finally:
            self.root.after(0, lambda: self.send_btn.config(state=tk.NORMAL))

    def recover_session(self):
        data = self.storage.load()
        if data:
            self.history = data
            for msg in self.history:
                tag = "user" if msg['role'] == 'user' else None
                label = "Tú: " if msg['role'] == 'user' else "Habi: "
                self.append_text(f"{label}{msg['content']}\n\n", tag)

if __name__ == "__main__":
    app_root = tk.Tk()
    HabiApp(app_root, OllamaClient(), SessionManager())
    app_root.mainloop()
