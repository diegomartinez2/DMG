#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  untitled.py
#
#  Copyright 2024 Diego Martinez Gutierrez <diego.martinez@ehu.eus>
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
import pygame
pygame.init()

width = 800
height = 600
screen = pygame.display.set_mode((width, height))

white = (255, 255, 255)
red = (255, 0, 0)
ball_pos = [width / 2, height / 2]
ball_radius = 20
speed = [2, 2]
gravity = 0.5
dragging = False
stored_speed = None

running = True
while running:
    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            running = False
        elif event.type == pygame.MOUSEBUTTONDOWN:
            x, y = event.pos
            dx = x - ball_pos[0]
            dy = y - ball_pos[1]
            distance = dx**2 + dy**2
            if distance <= ball_radius**2:
                dragging = True
                stored_speed = list(speed) # Store the current speed
        elif event.type == pygame.MOUSEBUTTONUP:
            dragging = False
            speed = stored_speed # Restore the stored speed

    if not dragging:
        # Apply gravity
        speed[1] -= gravity

        ball_pos[0] += speed[0]
        ball_pos[1] += speed[1]

        # Bounce off the sides of the screen
        if ball_pos[0] + ball_radius > width or ball_pos[0] - ball_radius < 0:
            speed[0] *= -1

        # Bounce off the top of the screen
        if ball_pos[1] - ball_radius < 0:
            speed[1] *= -1

        # Bounce off the bottom of the screen
        if ball_pos[1] + ball_radius > height:
            ball_pos[1] = height - ball_radius
            speed[1] *= -1

    if dragging:
        x, y = pygame.mouse.get_pos()
        ball_pos = [x, y]

    screen.fill(white)
    pygame.draw.circle(screen, red, (int(ball_pos[0]), int(ball_pos[1])), ball_radius)
    pygame.display.flip()

pygame.quit()
