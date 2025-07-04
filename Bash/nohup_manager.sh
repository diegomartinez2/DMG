#!/bin/bash
# Este script de bash busca la ejecucion de nohup mediante 'ps aux' y permite
# pausar el calculo o continuarlo.
# NOTA: no funciona como esperaba ya que los programas lanzados con nohup no
# se muestran como 'nohup' el el 'ps aux'
# --- Colores para la salida en terminal ---
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# --- Función para mostrar mensajes de error ---
function error_message() {
    echo -e "${RED}ERROR: $1${NC}" >&2
}

# --- Función para mostrar mensajes de éxito ---
function success_message() {
    echo -e "${GREEN}ÉXITO: $1${NC}"
}

# --- Función para listar los procesos nohup ---
function list_nohup_processes() {
    echo -e "${BLUE}--- Procesos 'nohup' en ejecución ---${NC}"
    # Usamos ps aux para obtener todos los procesos
    # filtramos por 'nohup ' (con espacio al final para evitar falsos positivos)
    # filtamos 'grep' de la lista de resultados
    # awk para formatear la salida: PID y comando
    # sed para eliminar el 'nohup ' del comando mostrado
    # grep -v "^$" para eliminar líneas vacías
    # sort -nk1 para ordenar por PID numéricamente
    # Esto crea un array asociativo (diccionario) 'processes' donde la clave es el número de la lista y el valor es el PID.
    # También imprime la lista formateada.

    declare -A processes # Array asociativo para mapear el número de lista al PID
    local i=1
    # Captura la salida de ps y la procesa línea por línea
    ps aux | grep "nohup " | grep -v "grep" | grep -v "$0" | awk '{print $2, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20}' | while read -r pid cmd; do
        # Elimina 'nohup ' del inicio del comando para una mejor visualización
        clean_cmd=$(echo "$cmd" | sed 's/^nohup //')
        # Si el comando es demasiado largo, lo truncamos para que la tabla sea legible
        if [[ ${#clean_cmd} -gt 70 ]]; then
            clean_cmd="${clean_cmd:0:67}..."
        fi
        printf "  ${YELLOW}%-3s${NC} | PID: ${BLUE}%-7s${NC} | Comando: ${GREEN}%s${NC}\n" "$i" "$pid" "$clean_cmd"
        processes[$i]="$pid"
        ((i++))
    done

    # Exporta el array para que sea accesible fuera de la función (en el ámbito del script)
    # Esto es crucial para que el bucle principal pueda acceder a los PIDs seleccionados
    # Si bash es versión 4 o superior, se puede usar declare -gA
    # Para compatibilidad, lo pasaremos como una variable global de esta forma.
    _processes_pids_map=$(declare -p processes)
}

# --- Función principal del menú ---
function main_menu() {
    local choice
    local selected_pid
    local action

    while true; do
        list_nohup_processes
        # Evalúa el string del array asociativo para restaurarlo en el ámbito actual
        eval "$_processes_pids_map"

        if [[ ${#processes[@]} -eq 0 ]]; then
            echo -e "${YELLOW}No se encontraron procesos 'nohup' en ejecución.${NC}"
        fi

        echo -e "\n${BLUE}--- Opciones ---${NC}"
        echo "1. Pausar un proceso (SIGSTOP)"
        echo "2. Continuar un proceso (SIGCONT)"
        echo "3. Recargar lista de procesos"
        echo "4. Salir"
        echo -e "${BLUE}----------------${NC}"
        read -p "Elige una opción (1-4): " choice

        case $choice in
            1|2)
                read -p "Introduce el número del proceso a ${choice_text}: " process_num
                if [[ -z "${processes[$process_num]}" ]]; then
                    error_message "Número de proceso inválido o no encontrado."
                    continue
                fi
                selected_pid="${processes[$process_num]}"

                if [[ "$choice" -eq 1 ]]; then
                    action="SIGSTOP"
                    choice_text="pausar"
                else
                    action="SIGCONT"
                    choice_text="continuar"
                fi

                echo -e "Intentando ${choice_text} el proceso PID: ${BLUE}$selected_pid${NC}..."
                if kill -"$action" "$selected_pid" 2>/dev/null; then
                    success_message "Proceso PID $selected_pid ${choice_text}do."
                else
                    error_message "No se pudo ${choice_text} el proceso PID $selected_pid. ¿Permisos insuficientes o proceso ya terminado?"
                fi
                ;;
            3)
                echo -e "${YELLOW}Recargando lista de procesos...${NC}"
                ;;
            4)
                echo -e "${BLUE}Saliendo del gestor de procesos.${NC}"
                exit 0
                ;;
            *)
                error_message "Opción inválida. Por favor, elige un número entre 1 y 4."
                ;;
        esac
        echo "" # Línea en blanco para mejor legibilidad
        sleep 1 # Pequeña pausa para que el usuario vea el mensaje
    done
}

# --- Iniciar el menú principal ---
main_menu
