#!/bin/bash

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

# --- Función para listar los procesos 'chimes' del usuario actual ---
function list_chimes_processes() {
    echo -e "${BLUE}--- Procesos 'chimes' en ejecución (solo usuario actual) ---${NC}"
    # Usamos ps -u $(whoami) para obtener solo los procesos del usuario actual.
    # Filtramos por 'chimes' en el comando.
    # Filtramos 'grep' de la lista de resultados y el propio script.
    # awk para formatear la salida: PID y el comando completo.
    # sort -nk1 para ordenar por PID numéricamente.

    declare -A processes # Array asociativo para mapear el número de lista al PID
    local i=1
    # Captura la salida de ps y la procesa línea por línea
    ps -u $(whoami) | grep "chimes" | grep -v "grep" | grep -v "$0" | awk '{print $2, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20}' | while read -r pid cmd; do
        # Si el comando es demasiado largo, lo truncamos para que la tabla sea legible
        if [[ ${#cmd} -gt 70 ]]; then
            clean_cmd="${cmd:0:67}..."
        else
            clean_cmd="$cmd"
        fi
        printf "  ${YELLOW}%-3s${NC} | PID: ${BLUE}%-7s${NC} | Comando: ${GREEN}%s${NC}\n" "$i" "$pid" "$clean_cmd"
        processes[$i]="$pid"
        ((i++))
    done

    # Exporta el array para que sea accesible fuera de la función (en el ámbito del script)
    _processes_pids_map=$(declare -p processes)
}

# --- Función principal del menú ---
function main_menu() {
    local choice
    local selected_pid
    local action
    local choice_text # Declarar aquí para asegurar que esté definida

    while true; do
        list_chimes_processes
        # Evalúa el string del array asociativo para restaurarlo en el ámbito actual
        eval "$_processes_pids_map"

        if [[ ${#processes[@]} -eq 0 ]]; then
            echo -e "${YELLOW}No se encontraron procesos 'chimes' en ejecución para el usuario actual.${NC}"
        fi

        echo -e "\n${BLUE}--- Opciones ---${NC}"
        echo "1. Pausar un proceso (SIGSTOP)"
        echo "2. Continuar un proceso (SIGCONT)"
        echo "3. Recargar lista de procesos"
        echo "4. Salir"
        echo -e "${BLUE}----------------${NC}"
        read -p "Elige una opción (1-4): " choice

        case $choice in
            1) # Pausar
                choice_text="pausar"
                read -p "Introduce el número del proceso a ${choice_text}: " process_num
                if [[ -z "${processes[$process_num]}" ]]; then
                    error_message "Número de proceso inválido o no encontrado."
                    continue
                fi
                selected_pid="${processes[$process_num]}"
                action="SIGSTOP"
                echo -e "Intentando ${choice_text} el proceso PID: ${BLUE}$selected_pid${NC}..."
                if kill -"$action" "$selected_pid" 2>/dev/null; then
                    success_message "Proceso PID $selected_pid ${choice_text}do."
                else
                    error_message "No se pudo ${choice_text} el proceso PID $selected_pid. ¿Permisos insuficientes o proceso ya terminado?"
                fi
                ;;
            2) # Continuar
                choice_text="continuar"
                read -p "Introduce el número del proceso a ${choice_text}: " process_num
                if [[ -z "${processes[$process_num]}" ]]; then
                    error_message "Número de proceso inválido o no encontrado."
                    continue
                fi
                selected_pid="${processes[$process_num]}"
                action="SIGCONT"
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
