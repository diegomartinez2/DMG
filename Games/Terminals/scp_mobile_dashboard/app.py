#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  app.py
#
#  Copyright 2025 Diego Martinez Gutierrez <diego.martinez@ehu.eus>
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#
#
# app.py
# Fundación SCP – Dashboard Móvil 2025 (Flask + WebSockets)

from flask import Flask, render_template, request
from flask_socketio import SocketIO, emit
import threading
import time
import random

app = Flask(__name__)
app.config['SECRET_KEY'] = 'SCP-001-O5-7-CLASSIFIED'
socketio = SocketIO(app, cors_allowed_origins="*", async_mode='eventlet')

# Estado global
alert_level = 2
alerts = ["WHITE", "GREEN", "BLUE", "AMBER", "RED", "KETER"]
keter_mode = False

# Simulación de eventos anómalos
def anomaly_simulator():
    global alert_level, keter_mode
    while True:
        time.sleep(random.randint(15, 60))
        if keter_mode: continue

        event = random.choice([
            ("SCP-173", "Blink detection failed – possible breach"),
            ("SCP-096", "Image SCP-096-1 detected on camera 4"),
            ("SCP-079", "AI attempted network escape"),
            ("SCP-682", "Containment acid neutralized – adaptation complete"),
        ])
        socketio.emit('anomaly', {
            'scp': event[0],
            'message': event[1],
            'level': 'CRITICAL'
        })

        if random.random() < 0.2:
            alert_level = min(5, alert_level + 1)
            socketio.emit('alert_change', {'level': alert_level, 'name': alerts[alert_level]})

threading.Thread(target=anomaly_simulator, daemon=True).start()

@app.route('/')
def index():
    return render_template('index.html')

@socketio.on('connect')
def handle_connect():
    emit('status', {
        'alert': alerts[alert_level],
        'level': alert_level,
        'site': 'Site-19',
        'keter': keter_mode
    })

@socketio.on('execute_scp')
def handle_scp(data):
    scp = data['scp']
    if scp == "SCP-001":
        global keter_mode, alert_level
        keter_mode = True
        alert_level = 5
        emit('keter_activated', broadcast=True)
        emit('anomaly', {
            'scp': 'SCP-001',
            'message': 'XK-CLASS END-OF-WORLD SCENARIO INITIATED',
            'level': 'CRITICAL'
        }, broadcast=True)
    else:
        emit('log', {
            'scp': scp,
            'message': f"{scp} containment protocols verified",
            'level': 'SUCCESS'
        })

if __name__ == '__main__':
    print("SCP Mobile Dashboard 2025 – ONLINE")
    print("Accede desde tu móvil: http://TU_IP:5000")
    socketio.run(app, host='0.0.0.0', port=5000)
