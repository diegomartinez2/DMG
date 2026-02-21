# Control Panel Simulation - Applying SOLID Principles

import time  # For delays
import random # For simulating probability
import tkinter as tk

# ------------------ Component Classes ------------------

class Button:
    """
    Represents a button with states (off, on, testing) and a callback function.
    Implements SRP: Single Responsibility - Manages button state and UI update.
    """
    def __init__(self, label, callback):
        self.label = label
        self.state = "off"  # Initial state
        self.callback = callback  # Function to call when state changes
        self.color = "white" #Initial color

    def press(self):
        """Simulates button press, updates state and calls callback."""
        self.state = "on"
        self.color = "white"
        self.callback() # Call the associated callback function


class Switch:
    """Represents a multi-phase switch with configurable states and behaviors."""
    def __init__(self, phase, callback):
        self.phase = phase
        self.state = "off"  # Initial state
        self.callback = callback

    def operate(self):
        """Simulates operating the switch, updates state and calls callback."""
        self.state = "on"
        self.callback()


class PlasmaControl:
    """Manages Plasma injection and purging simulation."""
    def __init__(self):
        self.injection_running = False
        self.purge_running = False

    def start_injection(self):
        """Starts plasma injection simulation with a 50% chance of error."""
        if not self.injection_running:
            self.injection_running = True
            print("Plasma injection started...")
            time.sleep(10)  # Simulate 10 second duration
            error_chance = 0.5
            if random.random() < error_chance:
                print("Plasma injection failed! (Red light)")
                self.set_button_color("red")
            else:
                print("Plasma injection successful.")
                self.set_button_color("green")

    def stop_injection(self):
        self.injection_running = False

    def start_purge(self):
        """Starts plasma purge simulation with a 50% chance of error."""
        if not self.purge_running:
            self.purge_running = True
            print("Plasma purge started...")
            time.sleep(10)  # Simulate 10 second duration
            error_chance = 0.5
            if random.random() < error_chance:
                print("Plasma purge failed! (Red light)")
                self.set_button_color("red")
            else:
                print("Plasma purge successful.")
                self.set_button_color("green")
        self.stop_purge()

    def stop_purge(self):
        self.purge_running = False

    def set_button_color(self, color):
        """Sets the color of the buttons."""
        print(f"Button color set to: {color}")  #Simulate updating button color

def create_gui():
    global window, button1, button2, button3, plasma_control

    window = tk.Tk()
    window.title("Plasma Control Panel")

    #Button 1: Injection
    button1 = tk.Button(window, text="Injection", width=10, height=2, command=lambda: plasma_control.start_injection())
    button1.grid(row=0, column=0)

    #Button 2: Purge
    button2 = tk.Button(window, text="Purge", width=10, height=2, command=lambda: plasma_control.start_purge())
    button2.grid(row=0, column=1)

    #Button 3: Testing
    button3 = tk.Button(window, text="Testing", width=10, height=2, command=lambda: print("Testing function called"))
    button3.grid(row=1, column=0)

    #Plasma Control (this remains the same)
    plasma_control = PlasmaControl()

    window.mainloop()


# ------------------ Example Usage ------------------
def button_press_callback():
    print("Button pressed!")
    #Simulate a change in the button state
    pass

# Create Button instances
button1 = Button("Injection", button_press_callback)
button2 = Button("Purge", button_press_callback)
button3 = Button("Testing", button_press_callback)

#Create Plasma Control
plasma_control = PlasmaControl()


# Example scenario - simulating button presses and plasma operations
plasma_control.start_injection()
button1.press() #Simulate a button press
plasma_control.start_purge()
button2.press() #Simulate a button press
