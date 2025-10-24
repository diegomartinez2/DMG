from ursina import *

app = Ursina()

# Configura el escenario
escenario = Entity(
    model='escenario.obj',  # Carga el modelo del terreno
    texture='escenario_textura.png',  # Textura del escenario
    scale=(1, 1, 1),  # Ajusta según el tamaño del modelo
    collider='mesh'  # Colisionador basado en el modelo
)

# Configura el jugador
jugador = Entity(
    model='jugador.obj',  # Carga el modelo del jugador
    texture='jugador_textura.png',  # Textura del jugador
    scale=0.5,  # Escala pequeña para el personaje
    position=(0, 2, 0),  # Posición inicial sobre el escenario
    collider='box'  # Colisionador simple para el jugador
)

# Variables para el movimiento y salto
velocidad = 5  # Velocidad de movimiento
velocidad_salto = 5  # Fuerza del salto
gravedad = -20  # Gravedad para el salto
velocidad_y = 0  # Velocidad vertical (para salto)
en_suelo = True  # Estado para verificar si está en el suelo

# Configura la cámara en tercera persona
camara = Entity(parent=jugador, position=(0, 2, -5))  # Cámara detrás y arriba del jugador
camera.look_at(jugador)  # La cámara siempre mira al jugador

# Función para ajustar la altura al terreno
def ajustar_altura_al_terreno():
    global en_suelo, velocidad_y
    # Lanza un rayo desde el jugador hacia abajo
    rayo = raycast(jugador.position, direction=(0, -1, 0), distance=10, ignore=[jugador])
    if rayo.hit and rayo.entity == escenario:
        altura_terreno = rayo.world_point.y
        # Si el jugador está cerca del terreno, ajústalo
        if abs(jugador.y - altura_terreno) < 0.5:
            jugador.y = altura_terreno + 0.5  # Ajusta altura (0.5 es el offset del modelo)
            en_suelo = True
            velocidad_y = 0
        elif jugador.y > altura_terreno + 0.5:
            en_suelo = False
    else:
        en_suelo = False

# Función de actualización para movimiento y física
def update():
    global velocidad_y, en_suelo

    # Movimiento horizontal con WASD
    direccion = Vec3(
        (held_keys['d'] - held_keys['a']) * velocidad,
        0,
        (held_keys['w'] - held_keys['s']) * velocidad
    ).normalized() * velocidad * time.dt
    jugador.position += direccion
    if direccion.length() > 0:
    jugador.rotation_y = math.atan2(direccion.x, direccion.z) * 180 / math.pi

    # Aplicar gravedad
    if not en_suelo:
        velocidad_y += gravedad * time.dt
        jugador.y += velocidad_y * time.dt

    # Ajustar la altura al terreno
    ajustar_altura_al_terreno()

    # Actualizar posición de la cámara (permanece fija relativa al jugador)
    camera.world_position = jugador.world_position + Vec3(0, 2, -5)
    camera.look_at(jugador)
    camera.rotation_y += mouse.velocity[0] * 100
    camera.rotation_x -= mouse.velocity[1] * 100

# Función para manejar la entrada (salto)
def input(key):
    global velocidad_y, en_suelo
    if key == 'space' and en_suelo:
        velocidad_y = velocidad_salto
        en_suelo = False

# Habilitar una cámara de editor para pruebas (opcional, desactívala si quieres)
# EditorCamera()

# Iniciar el juego
app.run()
