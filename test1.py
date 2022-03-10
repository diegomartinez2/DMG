# pyODE example 2: Connecting bodies with joints

import pygame
from pygame.locals import *
import ode
import numpy


def coord(x,y):
    "Convert world coordinates to pixel coordinates."
    return 320+170*x, 400-170*y


# Initialize pygame
pygame.init()

t=0.0

# Open a display
srf = pygame.display.set_mode((640,480))

# Simulation loop...

fps = 50
dt = 1.0/fps
loopFlag = True
clk = pygame.time.Clock()

while loopFlag:
    events = pygame.event.get()
    for e in events:
        if e.type==QUIT:
            loopFlag=False
        if e.type==KEYDOWN:
            loopFlag=False

    #movimiento de la bola en la pantalla
    x=numpy.sin(t)*100+300
    y=numpy.cos(t)*100+200
    # Clear the screen
    srf.fill((255,255,255))

    # Draw the two bodies
    pygame.draw.circle(srf, (55,0,200), (x,y), 20, 0)

    pygame.display.flip()

    # Next simulation step
    t=t+dt
    # Try to keep the specified framerate
    clk.tick(fps)
