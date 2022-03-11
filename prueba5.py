# no funciona
import pygame
from pygame.locals import *
import numpy

# Initialize pygame
pygame.init()

# inicializa el tiempo
t=0.0

# Open a display
srf = pygame.display.set_mode((640,480))

# Simulation loop...

fps = 50
dt = 1.0/fps
g=10
x1=200
y1=200
x2=100
y2=200
vx1=0
vy1=-0.1
vx2=0
vy2=15
m1=1000
m2=10
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
    # r(t+h)=r(t)+h*v(t)+h**2*(F(t)/m)/2
    # v(t+h)=v(t)+h*(F(t)/m+F(t+h)/m)/2
    for indice in [1,2,3,4,5,6,7,8,9,10]:
	d=(x1-x2)**2+(y1-y2)**2
	x1=x1+(dt/10)*vx1+((dt/10)**2)*(g*m2/d)/2 *numpy.sign(x2-x1)
	y1=y1+(dt/10)*vy1+((dt/10)**2)*(g*m2/d)/2 *numpy.sign(y2-y1)
	x2=x2+(dt/10)*vx2+((dt/10)**2)*(g*m1/d)/2 *numpy.sign(x1-x2)
	y2=y2+(dt/10)*vy2+((dt/10)**2)*(g*m1/d)/2 *numpy.sign(y1-y2)
	dn=(x1-x2)**2+(y1-y2)**2
	vx1=vx1+(dt/10)*((g*m2/d)+(g*m2/dn))/2 *numpy.sign(x2-x1)    
	vy1=vy1+(dt/10)*((g*m2/d)+(g*m2/dn))/2 *numpy.sign(y2-y1)    
	vx2=vx2+(dt/10)*((g*m1/d)+(g*m1/dn))/2 *numpy.sign(x1-x2)    
	vy2=vy2+(dt/10)*((g*m1/d)+(g*m1/dn))/2 *numpy.sign(y1-y2)    

    # Clear the screen
    srf.fill((255,255,255))

    # Dibuja el circulo
    pygame.draw.circle(srf, (55,0,200), (x1,y1), numpy.sqrt(m1), 0)
    pygame.draw.circle(srf, (55,0,0), (x2,y2), numpy.sqrt(m2), 0)

    pygame.display.flip()

    # Next simulation step
    t=t+dt
    # Try to keep the specified framerate    
    clk.tick(fps)
