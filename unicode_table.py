import tkinter as tk
from tkinter import ttk
import unicodedata

def mostrar_unicode():
    root = tk.Tk()
    root.title("Navegador de Tabla Unicode y Atajos Ctrl+Shift+U")
    root.geometry("650x500")

    # Contenedor de búsqueda
    frame_busqueda = ttk.Frame(root, padding=10)
    frame_busqueda.pack(fill=tk.X)

    ttk.Label(frame_busqueda, text="Buscar carácter o nombre:").pack(side=tk.LEFT, padx=5)
    entry_busqueda = ttk.Entry(frame_busqueda)
    entry_busqueda.pack(side=tk.LEFT, fill=tk.X, expand=True, padx=5)

    # Tabla para mostrar los datos
    columns = ("char", "hex", "shortcut", "name")
    tree = ttk.Treeview(root, columns=columns, show="headings", selectmode="browse")
    
    tree.heading("char", text="Carácter")
    tree.heading("hex", text="Código Hex")
    tree.heading("shortcut", text="Secuencia (Ctrl+Shift+U)")
    tree.heading("name", text="Nombre Unicode")

    tree.column("char", width=80, anchor="center")
    tree.column("hex", width=100, anchor="center")
    tree.column("shortcut", width=180, anchor="center")
    tree.column("name", width=250, anchor="w")

    scrollbar = ttk.Scrollbar(root, orient=tk.VERTICAL, command=tree.yview)
    tree.configure(yscroll=scrollbar.set)
    
    # Corregido: padx y pady en lugar de padding
    tree.pack(side=tk.LEFT, fill=tk.BOTH, expand=True, padx=(10, 0), pady=10)
    scrollbar.pack(side=tk.RIGHT, fill=tk.Y, padx=(0, 10), pady=10)

    # Cargar un rango representativo de la tabla Unicode
    def cargar_caracteres(filtro=""):
        tree.delete(*tree.get_children())
        filtro = filtro.lower()
        
        rangos = [(32, 127), (160, 895), (8192, 11263), (128512, 128591)]
        
        for inicio, fin in rangos:
            for codepoint in range(inicio, fin):
                try:
                    char = chr(codepoint)
                    name = unicodedata.name(char, "DESCONOCIDO")
                    hex_code = f"{codepoint:04X}"
                    shortcut = f"Ctrl+Shift+U + {hex_code.lstrip('0')} + Enter"
                    
                    if filtro in char.lower() or filtro in name.lower() or filtro in hex_code.lower():
                        tree.insert("", tk.END, values=(char, f"U+{hex_code}", shortcut, name))
                except Exception:
                    continue

    def al_buscar(event):
        cargar_caracteres(entry_busqueda.get())

    entry_busqueda.bind("<KeyRelease>", al_buscar)

    # Carga inicial
    cargar_caracteres()

    root.mainloop()

if __name__ == "__main__":
    mostrar_unicode()
