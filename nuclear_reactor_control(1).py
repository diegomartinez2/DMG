import curses
import time
import random
import numpy as np

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

def draw_plot_window(win, voltage_data, time_data):
    win.clear()
    win.border()
    win.addstr(0, 2, "Voltage Plot (kV)", curses.A_BOLD)
    
    # Window dimensions for plotting
    height, width = win.getmaxyx()
    plot_height = height - 3  # Subtract border and title
    plot_width = width - 3   # Subtract borders
    max_voltage = 20.0  # Max voltage in kV
    
    # Clear plot if it reaches the end
    if len(voltage_data) >= plot_width:
        voltage_data.clear()
        time_data.clear()
    
    # Draw Y-axis (voltage scale)
    for i in range(plot_height):
        voltage = max_voltage * (plot_height - i - 1) / (plot_height - 1)
        if i % 4 == 0:  # Show labels every 4 lines
            win.addstr(i + 1, 1, f"{voltage:4.1f}")
    
    # Draw plot
    for x, voltage in enumerate(voltage_data):
        if x >= plot_width:
            break
        # Scale voltage to plot height
        y = int((voltage / max_voltage) * (plot_height - 1))
        if 0 <= y < plot_height:
            win.addstr(plot_height - y, x + 2, "*", curses.color_pair(1))
    
    # Draw X-axis label
    win.addstr(plot_height + 1, 2, "Time")
    
    win.refresh()

def update_reactor(reactor_on, temperature, power_level, turbine_rpm, watts_generated, voltage):
    current_time = time.time()
    
    if reactor_on:
        # Simulate temperature increase based on power level
        temperature += power_level * 0.1 + random.uniform(-5, 5)
        temperature -= 10  # Natural cooling
        temperature = max(20, min(temperature, 1000))
        power_level = max(0, min(power_level, 100))
        
        # Realistic turbine RPM model using logistic function
        max_rpm = 3600
        k = 0.01
        T_0 = 500
        turbine_rpm = max_rpm / (1 + np.exp(-k * (temperature - T_0)))
        
        watts_generated = (turbine_rpm / max_rpm) * 1000
        voltage = (watts_generated / 1000) * 20 * 0.95
    else:
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
    plot_win = stdscr.derwin(12, 40, 2, 46)

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
        draw_plot_window(plot_win, voltage_data, time_data)
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
                voltage_data.clear()
                time_data.clear()
        elif key == ord('q'):
            break

        temperature, power_level, turbine_rpm, watts_generated, voltage = update_reactor(
            reactor_on, temperature, power_level, turbine_rpm, watts_generated, voltage
        )

if __name__ == '__main__':
    curses.wrapper(main)