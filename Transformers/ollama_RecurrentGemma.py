import tkinter as tk
from tkinter import scrolledtext, messagebox
import threading
import ollama

class RecurrentGemmaApp:
    def __init__(self, root):
        self.root = root
        self.root.title("RecurrentGemma - Ollama GUI")
        self.root.geometry("700x600")
        self.root.configure(bg="#1e1e1e")

        # Configuración de estilos y fuentes
        self.font_main = ("Segoe UI", 11)
        self.bg_color = "#1e1e1e"
        self.text_color = "#e0e0e0"
        self.input_bg = "#2d2d2d"
        self.accent_color = "#4a90e2"

        # --- ÁREA SUPERIOR: SALIDA DE TEXTO ---
        self.output_label = tk.Label(root, text="Respuesta de RecurrentGemma:", bg=self.bg_color, fg=self.accent_color, font=("Segoe UI", 10, "bold"))
        self.output_label.pack(pady=(10, 0), padx=20, anchor="w")

        self.text_area = scrolledtext.ScrolledText(root, wrap=tk.WORD, font=self.font_main, bg=self.input_bg, fg=self.text_color, insertbackground="white", borderwidth=0)
        self.text_area.pack(padx=20, pady=10, fill=tk.BOTH, expand=True)
        self.text_area.config(state=tk.DISABLED) # Solo lectura inicialmente

        # --- ÁREA INFERIOR: ENTRADA DE PROMPT ---
        self.input_label = tk.Label(root, text="Introduce tu mensaje:", bg=self.bg_color, fg=self.accent_color, font=("Segoe UI", 10, "bold"))
        self.input_label.pack(pady=(10, 0), padx=20, anchor="w")

        self.input_area = tk.Text(root, height=4, font=self.font_main, bg=self.input_bg, fg=self.text_color, insertbackground="white", borderwidth=0)
        self.input_area.pack(padx=20, pady=10, fill=tk.X)
        self.input_area.bind("<Return>", self.handle_return) # Enviar con Enter

        # Botón de envío
        self.send_button = tk.Button(root, text="Enviar Prompt", command=self.send_prompt, bg=self.accent_color, fg="white", font=("Segoe UI", 10, "bold"), relief=tk.FLAT, cursor="hand2")
        self.send_button.pack(pady=(0, 20))

    def handle_return(self, event):
        # Evitar que el Enter cree una nueva línea y envíe el prompt
        if not event.state & 0x1: # Si no se pulsa Shift
            self.send_prompt()
            return "break"

    def update_output(self, text):
        """Añade texto al área de salida de forma segura desde hilos."""
        self.text_area.config(state=tk.NORMAL)
        self.text_area.insert(tk.END, text)
        self.text_area.see(tk.END)
        self.text_area.config(state=tk.DISABLED)

    def send_prompt(self):
        user_input = self.input_area.get("1.0", tk.END).strip()
        if not user_input:
            return

        # Limpiar entrada y mostrar el mensaje del usuario en la salida
        self.input_area.delete("1.0", tk.END)
        self.update_output(f"\nUsted: {user_input}\n\nAI: ")

        # Deshabilitar botón mientras procesa
        self.send_button.config(state=tk.DISABLED, text="Procesando...")

        # Ejecutar la llamada a Ollama en un hilo separado para no bloquear la UI
        thread = threading.Thread(target=self.call_ollama, args=(user_input,))
        thread.start()

    def call_ollama(self, prompt):
        try:
            # Usamos streaming para que la respuesta aparezca poco a poco
            response = ollama.chat(
                model='recurrent-gemma',
                messages=[{'role': 'user', 'content': prompt}],
                stream=True
            )

            for chunk in response:
                content = chunk['message']['content']
                # Actualizar la UI desde el hilo principal
                self.root.after(0, self.update_output, content)

            self.root.after(0, self.update_output, "\n" + "-"*30 + "\n")

        except Exception as e:
            error_msg = f"\nError al conectar con Ollama: {str(e)}\n"
            self.root.after(0, self.update_output, error_msg)
            self.root.after(0, lambda: messagebox.showerror("Error", "Asegúrate de que Ollama esté corriendo y el modelo 'recurrent-gemma' esté descargado."))

        finally:
            self.root.after(0, lambda: self.send_button.config(state=tk.NORMAL, text="Enviar Prompt"))

if __name__ == "__main__":
    # Asegúrate de tener instalada la librería: pip install ollama
    root = tk.Tk()
    app = RecurrentGemmaApp(root)
    root.mainloop()
