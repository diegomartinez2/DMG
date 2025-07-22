#!/usr/bin/env python
"""
How to Use

    Navigation: Use the UP and DOWN arrow keys to select a button (Toggle Reactor, Increase Power, Decrease Power, Emergency Shutdown).
    Activate: Press ENTER to execute the selected action.
    Quit: Press 'q' to exit the program.
    Indicators:
        Reactor Status: Shows ON (green) or OFF (red).
        Temperature: Displays current temperature with a bar (turns yellow for warnings above 800°C, red for critical above 900°C).
        Power Level: Shows the power output percentage with a bar.
    Controls:
        Toggle Reactor: Turns the reactor ON or OFF.
        Increase/Decrease Power: Adjusts power level by 10% (only when reactor is ON).
        Emergency Shutdown: Turns off the reactor, resets power to 0%, and cools it to 20°C.

Simulation Details

    When the reactor is ON, temperature increases based on power level with some random fluctuation, and cools slightly each cycle.
    When OFF, the reactor cools down to a minimum of 20°C, and power drops to 0%.
    The interface updates every 100ms, creating a dynamic simulation.
    Alerts appear for high (yellow) or critical (red) temperatures.
"""
import curses
import time
import random

def draw_panel(stdscr, reactor_on, temperature, power_level, selected_button):
    stdscr.clear()
    height, width = stdscr.getmaxyx()

    # Ensure the terminal is large enough
    if height < 20 or width < 60:
        stdscr.addstr(0, 0, "Terminal too small! Resize to at least 60x20.")
        stdscr.refresh()
        return

    # Title
    title = "Nuclear Reactor Control Panel"
    stdscr.addstr(1, (width - len(title)) // 2, title, curses.A_BOLD)

    # Status indicators
    status = "ON" if reactor_on else "OFF"
    status_color = curses.color_pair(1) if reactor_on else curses.color_pair(2)
    stdscr.addstr(4, 10, f"Reactor Status: {status}", status_color)

    # Temperature gauge
    temp_str = f"Temperature: {temperature:.1f}°C"
    temp_color = curses.color_pair(1) if temperature < 800 else curses.color_pair(3)
    stdscr.addstr(6, 10, temp_str, temp_color)
    temp_bar = int(temperature / 1000 * 20)
    stdscr.addstr(7, 10, "[" + "#" * temp_bar + "-" * (20 - temp_bar) + "]")

    # Power level gauge
    power_str = f"Power Level: {power_level:.1f}%"
    stdscr.addstr(9, 10, power_str)
    power_bar = int(power_level / 100 * 20)
    stdscr.addstr(10, 10, "[" + "#" * power_bar + "-" * (20 - power_bar) + "]")

    # Alerts
    if temperature > 900:
        stdscr.addstr(12, 10, "ALERT: CRITICAL TEMPERATURE!", curses.color_pair(3) | curses.A_BOLD)
    elif temperature > 800:
        stdscr.addstr(12, 10, "WARNING: High Temperature", curses.color_pair(3))

    # Control buttons
    buttons = [
        ("Toggle Reactor", 0),
        ("Increase Power", 1),
        ("Decrease Power", 2),
        ("Emergency Shutdown", 3)
    ]
    for i, (label, btn_id) in enumerate(buttons):
        attr = curses.A_REVERSE if selected_button == btn_id else 0
        stdscr.addstr(14 + i, 10, f"[ {label} ]", attr)

    # Instructions
    stdscr.addstr(height - 2, 10, "Use UP/DOWN to select, ENTER to activate, 'q' to quit")
    stdscr.refresh()

def update_reactor(reactor_on, temperature, power_level):
    if reactor_on:
        # Simulate temperature increase based on power level
        temperature += power_level * 0.1 + random.uniform(-5, 5)
        # Natural cooling
        temperature -= 10
        # Ensure temperature stays within bounds
        temperature = max(20, min(temperature, 1000))
        # Power level stays within 0-100
        power_level = max(0, min(power_level, 100))
    else:
        # Cool down when off
        temperature -= 20
        temperature = max(20, temperature)
        power_level = 0
    return temperature, power_level

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
    selected_button = 0

    # Set up curses
    stdscr.timeout(100)  # Refresh every 100ms
    curses.curs_set(0)   # Hide cursor

    while True:
        draw_panel(stdscr, reactor_on, temperature, power_level, selected_button)

        try:
            key = stdscr.getch()
        except:
            key = -1

        # Handle input
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
        elif key == ord('q'):
            break

        # Update reactor state
        temperature, power_level = update_reactor(reactor_on, temperature, power_level)

if __name__ == '__main__':
    curses.wrapper(main)
