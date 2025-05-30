import curses
import time
import random
import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime
import os

# Global list to store voltage data for plotting
voltage_data = []
time_data = []

def draw_reactor_window(win, reactor_on, temperature, power_level, selected_button):
    win.clear()
    win.border()
    win.addstr(0, 2, "Reactor Control", curses.A_BOLD)
    
    # Status indicators
    status = "ON" if reactor_on else "OFF"
    status_color = curses.color_pair(1) if reactor_on else curses.color_pair(2)
    win.addstr(2, 2, f"Status: {status}", status_color)

    # Temperature gauge
    temp_str = f"Temperature: {temperature:.1f}°C"
    temp_color = curses.color_pair(1) if temperature < 800 else curses.color_pair(3)
    win.addstr(4, 2, temp_str, temp_color)
    temp_bar = int(temperature / 1000 * 15)
    win.addstr(5, 2, "[" + "#" * temp_bar + "-" * (15 - temp_bar) + "]")

    # Power level gauge
    power_str = f"Power Level: {power_level:.1f}%"
    win.addstr(7, 2, power_str)
    power_bar = int(power_level / 100 * 15)
    win.addstr(8, 2, "[" + "#" * power_bar + "-" * (15 - power_bar) + "]")

    # Alerts
    if temperature > 900:
        win.addstr(10, 2, "CRITICAL TEMP!", curses.color_pair(3) | curses.A_BOLD)
    elif temperature > 800:
        win.addstr(10, 2, "High Temp", curses.color_pair(3))

    win.refresh()

def draw_turbine_window(win, turbine_rpm, watts_generated, voltage):
    win.clear()
    win.border()
    win.addstr(0, 2, "Turbine Status", curses.A_BOLD)
    
    # Turbine indicators
    win.addstr(2, 2, f"RPM: {turbine_rpm:.0f}")
    rpm_bar = int(turbine_rpm / 3600 * 15)
    win.addstr(3, 2, "[" + "#" * rpm_bar + "-" * (15 - rpm_bar) + "]")

    win.addstr(5, 2, f"Power: {watts_generated:.0f} MW")
    watts_bar = int(watts_generated / 1000 * 15)
    win.addstr(6, 2, "[" + "#" * watts_bar + "-" * (15 - watts_bar) + "]")

    win.addstr(8, 2, f"Voltage: {voltage:.1f} kV")
    voltage_bar = int(voltage / 20 * 15)
    win.addstr(9, 2, "[" + "#" * voltage_bar + "-" * (15 - voltage_bar) + "]")

    win.refresh()

def draw_control_window(win, selected_button):
    win.clear()
    win.border()
    win.addstr(0, 2, "Controls", curses.A_BOLD)
    
    buttons = [
        ("Toggle Reactor", 0),
        ("Increase Power", 1),
        ("Decrease Power", 2),
        ("Emergency Shutdown", 3)
    ]
    for i, (label, btn_id) in enumerate(buttons):
        attr = curses.A_REVERSE if selected_button == btn_id else 0
        win.addstr(2 + i, 2, f"[ {label} ]", attr)
    
    win.addstr(7, 2, "UP/DOWN: Select, ENTER: Activate, q: Quit")
    win.refresh()

def update_plot(voltage_data, time_data):
    plt.figure(figsize=(8, 4))
    plt.plot(time_data, voltage_data, 'b-', label='Output Voltage (kV)')
    plt.title('Reactor Voltage Output Over Time')
    plt.xlabel('Time (s)')
    plt.ylabel('Voltage (kV)')
    plt.grid(True)
    plt.legend()
    plt.savefig('voltage_plot.png')
    plt.close()

def update_reactor(reactor_on, temperature, power_level, turbine_rpm, watts_generated, voltage):
    current_time = time.time()
    
    if reactor_on:
        # Simulate temperature increase based on power level
        temperature += power_level * 0.1 + random.uniform(-5, 5)
        temperature -= 10  # Natural cooling
        temperature = max(20, min(temperature, 1000))
        power_level = max(0, min(power_level, 100))
        
        # More realistic turbine RPM model based on temperature
        # Using a logistic function to model RPM with temperature
        max_rpm = 3600
        k = 0.01  # Steepness of the curve
        T_0 = 500  # Midpoint temperature for RPM scaling
        turbine_rpm = max_rpm / (1 + np.exp(-k * (temperature - T_0)))
        
        # Power output proportional to turbine RPM
        watts_generated = (turbine_rpm / max_rpm) * 1000
        # Voltage proportional to power output with some efficiency factor
        voltage = (watts_generated / 1000) * 20 * 0.95  # 95% efficiency
    else:
        # Cool down and reset turbine when off
        temperature -= 20
        temperature = max(20, temperature)
        power_level = 0
        turbine_rpm = 0
        watts_generated = 0
        voltage = 0

    # Store voltage data for plotting
    if len(voltage_data) == 0:
        start_time = current_time
    voltage_data.append(voltage)
    time_data.append(current_time - start_time)
    
    # Keep only last 100 points to prevent memory issues
    if len(voltage_data) > 100:
        voltage_data.pop(0)
        time_data.pop(0)
    
    # Update plot
    update_plot(voltage_data, time_data)
    
    return temperature, power_level, turbine_rpm, watts_generated, voltage

def main(stdscr):
    # Initialize colors
    curses.start_color()
    curses.init_pair(1, curses.COLOR_GREEN, curses.COLOR_BLACK)
    curses.init_pair(2, curses.COLOR_RED, curses.COLOR_BLACK)
    curses.init_pair(3, curses.COLOR_YELLOW, curses.COLOR_BLACK)

    # Initial reactor state
    reactor_on = False
    temperature = 20.0
    power_level = 0.0
    turbine_rpm = 0.0
    watts_generated = 0.0
    voltage = 0.0
    selected_button = 0

    # Create windows
    height, width = stdscr.getmaxyx()
    if height < 24 or width < 60:
        stdscr.addstr(0, 0, "Terminal too small! Resize to at least 60x24.")
        stdscr.refresh()
        time.sleep(1)
        return

    reactor_win = stdscr.derwin(12, 20, 2, 2)
    turbine_win = stdscr.derwin(12, 20, 2, 24)
    control_win = stdscr.derwin(10, 40, 14, 2)

    # Set up curses
    stdscr.timeout(100)
    curses.curs_set(0)

    # Title
    stdscr.addstr(0, (width - len("Nuclear Reactor Control Panel")) // 2, 
                 "Nuclear Reactor Control Panel", curses.A_BOLD)
    
    while True:
        draw_reactor_window(reactor_win, reactor_on, temperature, power_level, selected_button)
        draw_turbine_window(turbine_win, turbine_rpm, watts_generated, voltage)
        draw_control_window(control_win, selected_button)
        stdscr.refresh()

        try:
            key = stdscr.getch()
        except:
            key = -1

        if key == curses.KEY_UP:
            selected_button = (selected_button - 1) % 4
        elif key == curses.KEY_DOWN:
            selected_button = (selected_button + 1) % 4
        elif key == ord('\n'):
            if selected_button == 0:
                reactor_on = not reactor_on
            elif selected_button == 1 and reactor_on:
                power_level += 10
            elif selected_button == 2 and reactor_on:
                power_level -= 10
            elif selected_button == 3:
                reactor_on = False
                power_level = 0
                temperature = 20.0
                turbine_rpm = 0
                watts_generated = 0
                voltage = 0
        elif key == ord('q'):
            break

        temperature, power_level, turbine_rpm, watts_generated, voltage = update_reactor(
            reactor_on, temperature, power_level, turbine_rpm, watts_generated, voltage
        )

if __name__ == '__main__':
    curses.wrapper(main)