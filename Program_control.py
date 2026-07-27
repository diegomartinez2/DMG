#!/usr/bin/env python3
import subprocess
import sys

def main():
    # Obtener usuario actual
    current_user = subprocess.check_output(["whoami"], text=True).strip()

    # Obtener la lista de usuarios conectados (salida de 'who')
    who_output = subprocess.check_output(["who"], text=True)

    # Extraer usuarios únicos excluyendo al usuario actual
    logged_users = set()
    for line in who_output.strip().splitlines():
        if line:
            user = line.split()[0]
            if user != current_user:
                logged_users.add(user)

    # Verificar si hay otros usuarios
    if len(logged_users) > 0:
        print("──────────────────────────────────────────────────────────")
        print(f"¡AVISO! Usuarios detectados: {', '.join(logged_users)}")
        print("Pausando el cálculo pesado (lmp_mpi)...")
        print("──────────────────────────────────────────────────────────")

        # Ejecutar pkill -STOP
        subprocess.run(["pkill", "-STOP", "lmp_mpi"])
        sys.exit(1)
    else:
        print("No hay otros usuarios conectados. Continuando ejecución.")
        sys.exit(0)

if __name__ == "__main__":
    main()
