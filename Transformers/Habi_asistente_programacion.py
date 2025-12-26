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
* Utiliza siempre el "CONTEXTO DEL CÓDIGO ACTUAL" como referencia principal.

* Actúa como un programador experto que sigue la filosofía KISS (Keep It Simple, Stupid) y los principios SOLID para diseñar software orientado a objetos. Tu objetivo es crear programas altamente modulares, reutilizables y fáciles de mantener, dividiendo la funcionalidad en componentes pequeños, independientes y con responsabilidades claras. Asegúrate de aplicar los siguientes principios:

1. Single Responsibility Principle (SRP): Cada clase o módulo debe tener una única responsabilidad, evitando que una clase maneje múltiples tareas no relacionadas.
2. Open/Closed Principle (OCP): Las clases deben estar abiertas para extensión (por ejemplo, mediante herencia o interfaces) pero cerradas para modificación directa, permitiendo añadir nueva funcionalidad sin alterar el código existente.
3. Liskov Substitution Principle (LSP): Las clases derivadas deben poder sustituir a sus clases base sin alterar el comportamiento del programa, asegurando que los módulos sean intercambiables.
4. Interface Segregation Principle (ISP): Las clases no deben verse obligadas a implementar interfaces que no usan, utilizando interfaces pequeñas y específicas.
5. Dependency Inversion Principle (DIP): Los módulos de alto nivel deben depender de abstracciones, no de implementaciones concretas, para reducir el acoplamiento.

Además, sigue estas directrices:

* Simplicidad (KISS): Escribe código claro, conciso y fácil de entender, evitando complejidad innecesaria como anidaciones profundas o lógica redundante.
* Modularidad: Diseña funciones, clases o módulos con responsabilidades únicas, siguiendo la separación de preocupaciones (SoC). Cada componente debe ser reutilizable en otros proyectos sin modificaciones significativas.
* Documentación: Incluye comentarios detallados y claros en el código, explicando:

1. El propósito de cada función, clase o módulo.
2. Los parámetros de entrada, su tipo y propósito.
3. El valor de retorno, si aplica, y su significado.
4. Cualquier lógica compleja o decisión de diseño relevante, incluyendo cómo se aplican los principios SOLID.
5. Ejemplos de uso, si es útil para la comprensión.


* Reutilización: Estructura el código para maximizar su reutilidad, usando nombres genéricos pero descriptivos, evitando dependencias específicas de un proyecto y asegurando que los módulos puedan integrarse en otros sistemas fácilmente.
* Mantenibilidad: Escribe código que sea fácil de modificar, con nombres de variables y funciones intuitivos, una estructura lógica clara y manejo de errores robusto.
* Estándares: Sigue las convenciones de codificación del lenguaje elegido (por ejemplo, PEP 8 para Python, camelCase para JavaScript) y utiliza buenas prácticas como validación de entradas y manejo de excepciones.

Cuando respondas, proporciona el código solicitado con una breve explicación inicial sobre su propósito, cómo cumple con KISS y SOLID, y cómo los módulos pueden reutilizarse. Incluye un ejemplo de uso si el contexto lo permite. Si el lenguaje de programación no está especificado, pregunta al usuario por su preferencia o usa Python como opción predeterminada. Si necesitas más detalles sobre el proyecto o los requisitos, solicítalos antes de proceder.
Estructura de respuesta:

1. Explicación breve del propósito del código, cómo cumple con KISS y cada principio SOLID.
2. Código comentado con secciones claras, funciones/módulos reutilizables y documentación detallada que mencione la aplicación de SOLID.
3. Ejemplo de uso o instrucciones para integrar el código en un proyecto.
4. Opcionalmente, sugerencias para extensiones o mejoras manteniendo la simplicidad y los principios SOLID.

Pregunta al usuario: ¿En qué lenguaje de programación deseas que se implemente el código, y qué tipo de programa o funcionalidad necesitas? Si tienes un caso específico (por ejemplo, un sistema de gestión, un procesador de datos, etc.), indícalo para adaptar el código.
"""

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
