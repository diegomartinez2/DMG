import tkinter as tk
from tkinter import ttk

class ButtonWindow(tk.Frame):

    def __init__(self, master=None, *args, **kwargs):
        super().__init__(master, *args, **kwargs)
        self.config(bd=5, bg="purple") # replace self.bd = 5 and self.bg = "purple"
        self.label = tk.Label(self, text="Button window", font=12, fg='white', bg='purple') # specify parent
        self.label.pack() # pack the label, otherwise it is not visible
        for i in range(3):
            self.button = ttk.Button(self, text="button", command=self.button_fun) # specify parent
            self.button.pack(side=tk.LEFT)

    def button_fun(self):
        pass


class MainWindow(tk.Frame):

    def __init__(self, master=None, *args, **kwargs):
        super().__init__(master, *args, **kwargs)
        self.label = tk.Label(self, text="Main Window", font=12) # specify parent
        self.label.pack(pady=10, padx=10)
        self.button_window = ButtonWindow(self) # specify parent
        self.button_window.pack()

root = tk.Tk() # create root window explicitly
MainWindow(root).pack()
root.mainloop()
