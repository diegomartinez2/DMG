import pygame
import random
import math

pygame.init()
screen = pygame.display.set_mode((700, 700))

white = (255, 255, 255)
blue = (0, 0, 255)
yellow = (255, 255, 0)
grey = (200, 200, 200)
black = (0, 0, 0)

sun_radius = 50
center = (350, 350)

earth_x = 50
earth_y = 350

earth_orbit = 0
moon_orbit = 0

clock = pygame.time.Clock()

running = True

stars = [(random.randint(0, 699), random.randint(0, 699)) for x in range(140)]

while running:

    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            running = False

    # update earth position
    earth_x = math.cos(earth_orbit) * 300 + 350
    earth_y = -math.sin(earth_orbit) * 300 + 350

    # update moon position based on earth position
    moon_x = math.cos(moon_orbit) * 50 + earth_x
    moon_y = -math.sin(moon_orbit) * 50 + earth_y

    # update the moon and earth angles
    earth_orbit += .002
    moon_orbit += .01

    # reset the screen
    screen.fill(black)

    # draw the stars
    for star in stars:
        x, y = star[0], star[1]
        pygame.draw.line(screen, white, (x, y), (x, y))

    # draw the sun
    pygame.draw.circle(screen, yellow, center, sun_radius)

    # draw the earth
    pygame.draw.circle(screen, blue, (int(earth_x), int(earth_y)), 15)

    # draw the moon
    pygame.draw.circle(screen, grey, (int(moon_x), int(moon_y)), 5)

    pygame.display.flip()

    clock.tick(60)

pygame.quit()
