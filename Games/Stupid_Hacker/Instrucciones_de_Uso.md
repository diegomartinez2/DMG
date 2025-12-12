## **Guía para el Script Bash (hacker\_session.sh)**

Este script utiliza tmux (un multiplexor de terminales) para crear una sesión de terminal única dividida en 6 paneles, cada uno ejecutando un simulador de "hacking" con genact.

### **Requisitos**

1. **Sistema Operativo:** Linux o macOS (donde tmux y genact sean ejecutables).  
2. **tmux:** Debe estar instalado (sudo apt install tmux en Debian/Ubuntu).  
3. **genact:** El binario de genact debe estar instalado y accesible en tu $PATH.

### **Pasos para Ejecutar**

1. **Guarda el script:** Guarda el código proporcionado como hacker\_session.sh.  
2. **Da permisos de ejecución:** Abre tu terminal y ejecuta:  
   chmod \+x hacker\_session.sh

3. **Ejecuta el script:** Inicia la creación de la sesión:  
   ./hacker\_session.sh

   (El script imprimirá el nombre de la sesión: HackerSim).  
4. **Conéctate a la sesión:** Una vez que el script te lo indique, conéctate para ver la simulación:  
   tmux attach \-t HackerSim

### **Comandos Útiles dentro de Tmux**

* **Salir de la sesión (detach):** Ctrl+B y luego D. (La simulación sigue corriendo en segundo plano).  
* **Volver a la sesión:** tmux attach \-t HackerSim.  
* **Moverse entre paneles:** Ctrl+B y luego Flecha Arriba/Abajo/Izquierda/Derecha.  
* **Detener la sesión (cerrar todo):** Primero sal del programa (Ctrl+B y luego D), y luego ejecuta:  
  tmux kill-session \-t HackerSim

**Nota:** Si genact no está instalado, el script usará un comando while true básico para simular una salida de terminal en los paneles.
