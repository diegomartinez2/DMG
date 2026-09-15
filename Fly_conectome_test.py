from flygym import Fly, Camera
from flygym.simulation import Simulation

# 1. Instanciar la mosca biomecánica y el simulador de física
fly = Fly(enable_vision=True, enable_olfaction=True)
sim = Simulation(fly=fly, cameras=[Camera()])

# 2. Bucle de control sensoriomotor
obs, _ = sim.reset()
for _ in range(1000):
    # El modelo neuronal lee la visión y calcula los ángulos articulares de las patas
    action = mi_red_neuronal_conectoma.step(obs["vision"])
    obs, reward, terminated, truncated, info = sim.step(action)
