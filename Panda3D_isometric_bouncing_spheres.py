
from direct.showbase.ShowBase import ShowBase
from panda3d.core import NodePath, CollisionNode, CollisionSphere, CollisionPlane, Plane, Vec3
from panda3d.core import Point3, LVector3
from random import uniform
import sys

class IsometricGame(ShowBase):
    def __init__(self):
        ShowBase.__init__(self)

        # Deshabilitar control de cámara por ratón
        self.disableMouse()

        # Configurar cámara en vista isométrica
        self.camera.setPos(20, -20, 15)  # Posición para vista isométrica
        self.camera.lookAt(0, 0, 0)      # Apuntar al centro de la escena

        # Configurar luz ambiental para iluminar la escena
        ambient_light = self.render.attachNewNode("ambientLight")
        ambient_light.setLight(self.render, self.loader.loadModel("models/misc/light").node())
        self.render.setLight(ambient_light)

        # Crear caja delimitadora con planos de colisión
        self.setup_collision_walls()

        # Crear esferas
        self.spheres = []
        self.velocities = []
        self.num_spheres = 5
        self.setup_spheres()

        # Configurar sistema de colisiones
        self.setup_collisions()

        # Actualizar movimiento
        self.taskMgr.add(self.update, "update")

    def setup_collision_walls(self):
        # Definir límites de la caja (10x10x5)
        size = 5.0  # Media longitud de la caja
        planes = [
            Plane(Vec3(0, 0, 1), Point3(0, 0, -2.5)),  # Suelo (z = -2.5)
            Plane(Vec3(0, 0, -1), Point3(0, 0, 2.5)),  # Techo (z = 2.5)
            Plane(Vec3(1, 0, 0), Point3(-size, 0, 0)), # Pared izquierda
            Plane(Vec3(-1, 0, 0), Point3(size, 0, 0)), # Pared derecha
            Plane(Vec3(0, 1, 0), Point3(0, -size, 0)), # Pared trasera
            Plane(Vec3(0, -1, 0), Point3(0, size, 0)), # Pared frontal
        ]

        for plane in planes:
            collision_node = CollisionNode("wall")
            collision_node.addSolid(CollisionPlane(plane))
            self.render.attachNewNode(collision_node)

    def setup_spheres(self):
        # Cargar modelo de esfera
        sphere_model = self.loader.loadModel("models/misc/sphere")

        for i in range(self.num_spheres):
            # Crear esfera
            sphere = NodePath(f"sphere_{i}")
            sphere_model.copyTo(sphere)
            sphere.setScale(0.5)  # Tamaño de la esfera
            sphere.setPos(
                uniform(-4, 4),  # Posición inicial aleatoria en X
                uniform(-4, 4),  # Posición inicial aleatoria en Y
                uniform(-2, 2)   # Posición inicial aleatoria en Z
            )
            sphere.setColor(uniform(0, 1), uniform(0, 1), uniform(0, 1), 1)  # Color aleatorio
            sphere.reparentTo(self.render)

            # Añadir colisión a la esfera
            collision_node = CollisionNode(f"sphere_collider_{i}")
            collision_node.addSolid(CollisionSphere(0, 0, 0, 1))  # Radio ajustado a la escala
            sphere.attachNewNode(collision_node)

            # Velocidad inicial aleatoria
            velocity = LVector3(
                uniform(-2, 2),
                uniform(-2, 2),
                uniform(-1, 1)
            )
            self.spheres.append(sphere)
            self.velocities.append(velocity)

    def setup_collisions(self):
        # Configurar manejador de colisiones
        self.cTrav = self.collisionTraverser
        self.cHandler = self.collisionHandlerQueue
        for sphere in self.spheres:
            self.cTrav.addCollider(sphere.find("**/sphere_collider_*"), self.cHandler)

    def update(self, task):
        dt = globalClock.getDt()

        for i, sphere in enumerate(self.spheres):
            # Actualizar posición
            new_pos = sphere.getPos() + self.velocities[i] * dt
            sphere.setPos(new_pos)

            # Verificar colisiones
            for entry in self.cHandler.getEntries():
                from_node = entry.getFromNodePath()
                into_node = entry.getIntoNodePath()
                if from_node.getParent() == sphere:
                    # Colisión detectada
                    normal = entry.getSurfaceNormal(sphere)
                    self.velocities[i] = self.velocities[i] - 2 * self.velocities[i].dot(normal) * normal

        return task.cont

    def cleanup(self):
        self.taskMgr.remove("update")

# Iniciar el juego
app = IsometricGame()
app.run()

"""
Explicación del Código

    Clase Principal (IsometricGame):
        Hereda de ShowBase para inicializar Panda3D.
        Configura la cámara en una posición isométrica (20, -20, 15) mirando al origen.
    Caja Delimitadora:
        Se crean seis planos de colisión invisibles para formar una caja de 10x10x5 unidades.
        Los planos representan el suelo, el techo y las cuatro paredes laterales.
    Esferas:
        Se generan 5 esferas con posiciones iniciales aleatorias dentro de la caja.
        Cada esfera tiene un modelo cargado desde Panda3D, un color aleatorio y una esfera de colisión asociada.
        Las velocidades iniciales son aleatorias en los ejes X, Y y Z.
    Colisiones:
        Se usa el sistema de colisiones de Panda3D (CollisionTraverser y CollisionHandlerQueue).
        Cuando una esfera colisiona con una pared o otra esfera, se calcula la reflexión de su velocidad usando la normal de la superficie de colisión.
    Actualización:
        La función update mueve las esferas según sus velocidades y actualiza las trayectorias al detectar colisiones.
        Se usa globalClock.getDt() para mantener el movimiento independiente del framerate.

Cómo Ejecutar

    Requisitos:
        Instala Panda3D en Linux: pip install panda3d
        Asegúrate de tener Python 3.6+ instalado.
    Ejecutar:
        Guarda el script en un archivo, por ejemplo, isometric_bouncing_spheres.py.
        Ejecútalo con: python isometric_bouncing_spheres.py
    Resultado:
        Verás una ventana con una vista isométrica donde 5 esferas de colores rebotan dentro de una caja invisible, colisionando entre sí y con las paredes.

Notas

    Personalización: Puedes ajustar num_spheres, el tamaño de la caja (size), el tamaño de las esferas (setScale) o las velocidades iniciales (uniform) en el código.
    Limitaciones: El sistema de colisiones de Panda3D es básico. Para físicas más realistas (como fricción o elasticidad), podrías integrar Bullet Physics, pero esto añade complejidad.
    Mejoras posibles: Agregar texturas, sombras, controles de usuario (por ejemplo, mover la cámara con teclas) o un entorno visual para la caja.
    """
