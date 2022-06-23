import tkinter as tk

class ButtonWindow():

    def __init__(self, Frame, *args, **kwargs):
        self.label = tk.LabelFrame(Frame, text="Button window", font=12, bg = "purple")
        self.label.pack()
        for i in range(3):
            self.button = tk.Button(self.label, text="button", command= lambda: button_fun())
            self.button.pack(side = tk.LEFT)

def button_fun(self):
    pass


class MainWindow():

    def __init__(self, window, *args, **kwargs):

        myCoreFrame = tk.Frame(window)
        myCoreFrame.pack()

        self.label = tk.LabelFrame(myCoreFrame, text="Main Window", font=12, bg = "red")
        self.label.pack(pady=10,padx=10)
        self.button_window = ButtonWindow(self.label)

root = tk.Tk()
app = MainWindow(root)
root.mainloop()
