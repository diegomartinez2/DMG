"""
Explicación del código:

    Librerías usadas:
        Pygame: Para dibujar el reloj y manejar la ventana gráfica.
        datetime: Para obtener la hora actual del sistema.
        math: Para cálculos de ángulos (usando seno y coseno).
        asyncio y platform: Para compatibilidad con Pyodide en el navegador.
    Estructura:
        setup(): Configura la ventana del reloj.
        draw_clock(): Dibuja el círculo del reloj y las marcas de las horas.
        draw_hands(): Calcula los ángulos de las agujas (horas, minutos, segundos) según la hora actual y las dibuja.
        update_loop(): Actualiza la pantalla y maneja eventos (como cerrar la ventana).
        main(): Bucle principal asíncrono que controla la animación a 60 FPS.
    Detalles visuales:
        El reloj es un círculo con radio 150 píxeles.
        Las agujas tienen diferentes longitudes y grosores: horas (gruesa, corta), minutos (media), segundos (fina, larga, roja).
        Se dibujan 12 marcas para las horas, cada 30 grados.

Notas:

    El código está optimizado para ejecutarse en Pyodide, evitando operaciones de E/S locales y usando un bucle asíncrono.
    Si prefieres otra biblioteca o un enfoque diferente (por ejemplo, un reloj estático con Matplotlib o una interfaz web con HTML5 Canvas), por favor indícalos y puedo adaptar el código.
    Para probar el reloj, asegúrate de tener Pygame instalado (pip install pygame) si lo ejecutas localmente, o usa un entorno compatible con Pyodide para el navegador.
 """   
import pygame
import math
import datetime
import asyncio
import platform

# Inicialización de Pygame
pygame.init()
width, height = 400, 400
screen = pygame.display.set_mode((width, height))
clock = pygame.time.Clock()
FPS = 60

# Colores
WHITE = (255, 255, 255)
BLACK = (0, 0, 0)
RED = (255, 0, 0)
BLUE = (0, 0, 255)

# Configuración del reloj
center = (width // 2, height // 2)
radius = 150

def setup():
    pygame.display.set_caption("Reloj Analógico")

def draw_clock():
    # Fondo blanco
    screen.fill(WHITE)
    
    # Dibujar círculo del reloj
    pygame.draw.circle(screen, BLACK, center, radius, 2)
    
    # Dibujar marcas de hora
    for i in range(12):
        angle = math.radians(i * 30 - 90)  # 30 grados por hora
        x1 = center[0] + (radius - 10) * math.cos(angle)
        y1 = center[1] + (radius - 10) * math.sin(angle)
        x2 = center[0] + radius * math.cos(angle)
        y2 = center[1] + radius * math.sin(angle)
        pygame.draw.line(screen, BLACK, (x1, y1), (x2, y2), 2)

def draw_hands():
    # Obtener hora actual
    now = datetime.datetime.now()
    seconds = now.second
    minutes = now.minute + seconds / 60
    hours = now.hour % 12 + minutes / 60
    
    # Aguja de segundos (roja)
    second_angle = math.radians(seconds * 6 - 90)  # 6 grados por segundo
    second_length = radius - 20
    x = center[0] + second_length * math.cos(second_angle)
    y = center[1] + second_length * math.sin(second_angle)
    pygame.draw.line(screen, RED, center, (x, y), 2)
    
    # Aguja de minutos (negra)
    minute_angle = math.radians(minutes * 6 - 90)  # 6 grados por minuto
    minute_length = radius - 40
    x = center[0] + minute_length * math.cos(minute_angle)
    y = center[1] + minute_length * math.sin(minute_angle)
    pygame.draw.line(screen, BLACK, center, (x, y), 4)
    
    # Aguja de horas (negra)
    hour_angle = math.radians(hours * 30 - 90)  # 30 grados por hora
    hour_length = radius - 60
    x = center[0] + hour_length * math.cos(hour_angle)
    y = center[1] + hour_length * math.sin(hour_angle)
    pygame.draw.line(screen, BLACK, center, (x, y), 6)

def update_loop():
    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            return False
    
    draw_clock()
    draw_hands()
    pygame.display.flip()
    return True

async def main():
    setup()
    running = True
    while running:
        running = update_loop()
        await asyncio.sleep(1.0 / FPS)

if platform.system() == "Emscripten":
    asyncio.ensure_future(main())
else:
    if __name__ == "__main__":
        asyncio.run(main())
