import tkinter as tk
from tkinter import scrolledtext, messagebox, ttk
import threading
import ollama
import json
import os
from typing import List, Dict, Generator

# --- CONFIGURACIÓN ---
MODEL_NAME = "gemma3"

# Prompt de sistema optimizado para escritura creativa
HABI_SYSTEM_PROMPT = """Eres ‘Habi’, la musa de los escritores. 
Tu objetivo es ayudar a redactar novelas de alta calidad.
Utilizarás la información de contexto (Trama, Estilo, Personajes y Capítulos Anteriores) para que tus sugerencias sean coherentes.
* Si el usuario te pide 'escribir el capítulo', genera una narración rica y detallada.
* Si te pide un 'esquema', sé estructural.
* Mantén siempre el estilo solicitado por el autor.
* No utilices emojis. Responde en español."""

class SessionManager:
    def __init__(self, filename: str = "habi_novel_data.json"):
        self.filename = filename

    def save_all(self, history, sidebar_data):
        data = {
            "history": history,
            "sidebar": sidebar_data
        }
        try:
            with open(self.filename, 'w', encoding='utf-8') as f:
                json.dump(data, f, ensure_ascii=False, indent=4)
        except Exception as e:
            print(f"Error al guardar: {e}")

    def load_all(self):
        if os.path.exists(self.filename):
            try:
                with open(self.filename, 'r', encoding='utf-8') as f:
                    return json.load(f)
            except:
                return {"history": [], "sidebar": {}}
        return {"history": [], "sidebar": {}}

class OllamaClient:
    def __init__(self, model: str = MODEL_NAME):
        self.model = model

    def stream_chat(self, history: List[Dict[str, str]]) -> Generator[str, None, None]:
        try:
            response = ollama.chat(
                model=self.model,
                messages=history,
                stream=True
            )
            for chunk in response:
                yield chunk['message']['content']
        except Exception as e:
            raise ConnectionError(f"No se pudo conectar con Ollama: {e}")

class HabiApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Habi - Asistente de Escritura Creativa")
        self.root.geometry("1200x800")
        
        self.ai = OllamaClient()
        self.session = SessionManager()
        
        # Cargar datos previos
        saved_data = self.session.load_all()
        self.history = saved_data.get("history", [])
        if not self.history:
            self.history = [{"role": "system", "content": HABI_SYSTEM_PROMPT}]

        self.setup_ui()
        self.load_sidebar_content(saved_data.get("sidebar", {}))

    def setup_ui(self):
        # Contenedor Principal con divisor ajustable (Izquierda: Chat, Derecha: Herramientas)
        self.main_paned = tk.PanedWindow(self.root, orient=tk.HORIZONTAL, sashrelief=tk.RAISED, sashwidth=4)
        self.main_paned.pack(fill=tk.BOTH, expand=True)

        # --- LADO IZQUIERDO: CHAT ---
        self.chat_frame = tk.Frame(self.main_paned)
        self.main_paned.add(self.chat_frame, width=700)

        self.chat_display = scrolledtext.ScrolledText(self.chat_frame, wrap=tk.WORD, state=tk.DISABLED, bg="#fdfcfb", font=("Georgia", 11))
        self.chat_display.pack(padx=10, pady=10, fill=tk.BOTH, expand=True)
        self.chat_display.tag_configure("user", foreground="#2c3e50", font=("Georgia", 11, "bold"))

        self.input_frame = tk.Frame(self.chat_frame)
        self.input_frame.pack(fill=tk.X, padx=10, pady=(0, 10))

        self.input_text = tk.Text(self.input_frame, height=3, font=("Segoe UI", 10))
        self.input_text.pack(side=tk.LEFT, fill=tk.X, expand=True)
        self.input_text.bind("<Return>", lambda e: self.on_send() or "break")

        self.send_btn = tk.Button(self.input_frame, text="Preguntar a Habi", command=self.on_send, bg="#3498db", fg="white", padx=10)
        self.send_btn.pack(side=tk.RIGHT, padx=5, fill=tk.Y)

        self.status_label = tk.Label(self.chat_frame, text="Habi está lista para escribir contigo.", anchor=tk.W, fg="gray")
        self.status_label.pack(side=tk.BOTTOM, fill=tk.X, padx=10)

        # --- LADO DERECHO: CONTEXTO DE NOVELA ---
        self.sidebar_frame = tk.Frame(self.main_paned, bg="#f0f0f0")
        self.main_paned.add(self.sidebar_frame, width=450)

        # Usar un canvas con scroll para la barra lateral por si hay muchas ventanas
        self.setup_sidebar_widgets()

    def setup_sidebar_widgets(self):
        # Título General
        tk.Label(self.sidebar_frame, text="CONTEXTO DE LA NOVELA", font=("Segoe UI", 10, "bold"), bg="#f0f0f0").pack(pady=5)

        # 1. Trama del Capítulo
        self.create_context_box("Trama de este capítulo (Lo que debe pasar):", "trama_txt", height=8)
        
        # 2. Estilo de Escritura
        self.create_context_box("Estilo de Escritura (Ej: Primera persona, oscuro, poético):", "estilo_txt", height=3)

        # 3. Tarjetas de Personajes
        self.create_context_box("Tarjetas de Personajes (Nombres, rasgos, metas):", "personajes_txt", height=8)

        # 4. Capítulos Anteriores
        self.create_context_box("Capítulos Anteriores (Resumen o texto previo):", "anteriores_txt", height=12)

        # Botón para guardar manualmente
        tk.Button(self.sidebar_frame, text="Guardar todo el progreso", command=self.manual_save, bg="#27ae60", fg="white").pack(pady=10)

    def create_context_box(self, label_text, attr_name, height):
        tk.Label(self.sidebar_frame, text=label_text, bg="#f0f0f0", font=("Segoe UI", 9, "italic")).pack(anchor=tk.W, padx=10, pady=(10, 0))
        txt_area = scrolledtext.ScrolledText(self.sidebar_frame, height=height, font=("Consolas", 10), wrap=tk.WORD)
        txt_area.pack(fill=tk.X, padx=10, pady=2)
        setattr(self, attr_name, txt_area)

    def load_sidebar_content(self, sidebar_data):
        fields = {
            "trama": self.trama_txt,
            "estilo": self.estilo_txt,
            "personajes": self.personajes_txt,
            "anteriores": self.anteriores_txt
        }
        for key, widget in fields.items():
            content = sidebar_data.get(key, "")
            widget.insert("1.0", content)

    def get_sidebar_data(self):
        return {
            "trama": self.trama_txt.get("1.0", tk.END).strip(),
            "estilo": self.estilo_txt.get("1.0", tk.END).strip(),
            "personajes": self.personajes_txt.get("1.0", tk.END).strip(),
            "anteriores": self.anteriores_txt.get("1.0", tk.END).strip()
        }

    def append_text(self, content, tag=None):
        self.chat_display.config(state=tk.NORMAL)
        if tag:
            self.chat_display.insert(tk.END, content, tag)
        else:
            self.chat_display.insert(tk.END, content)
        self.chat_display.see(tk.END)
        self.chat_display.config(state=tk.DISABLED)

    def update_status(self, is_running):
        if is_running:
            self.status_label.config(text="Habi está escribiendo...", fg="#e67e22")
        else:
            self.status_label.config(text="Cambios guardados. Habi espera tus instrucciones.", fg="#27ae60")

    def manual_save(self):
        self.session.save_all(self.history, self.get_sidebar_data())
        messagebox.showinfo("Guardado", "Todo el progreso de la novela ha sido guardado.")

    def on_send(self):
        user_input = self.input_text.get("1.0", tk.END).strip()
        if not user_input: return

        # Recolectar contexto de las ventanas laterales
        context = self.get_sidebar_data()
        
        # Crear un mensaje "invisible" de contexto para la IA antes de su pregunta
        full_context_prompt = f"""
--- CONTEXTO ACTUAL DE LA NOVELA ---
TRAMA DEL CAPÍTULO: {context['trama']}
ESTILO: {context['estilo']}
PERSONAJES: {context['personajes']}
CAPÍTULOS ANTERIORES: {context['anteriores']}
-----------------------------------
INSTRUCCIÓN DEL USUARIO: {user_input}
"""
        self.input_text.delete("1.0", tk.END)
        
        # Añadimos a la historia visual el texto limpio del usuario, 
        # pero a la IA le enviamos el bloque con contexto.
        self.history.append({"role": "user", "content": full_context_prompt})
        self.append_text(f"Tú: {user_input}\n\nHabi: ", "user")

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
            
            # Guardar automáticamente después de cada respuesta
            self.session.save_all(self.history, self.get_sidebar_data())
            self.root.after(0, lambda: self.update_status(False))
            
        except Exception as e:
            self.root.after(0, lambda m=str(e): messagebox.showerror("Error", m))
        finally:
            self.root.after(0, lambda: self.send_btn.config(state=tk.NORMAL))

if __name__ == "__main__":
    root = tk.Tk()
    app = HabiApp(root)
    root.mainloop()
