def main():
    total_col1 = 0.0
    total_col2 = 0.0

    try:
        with open('DATOS.DAT', 'r') as file:
            for line in file:
                if line.strip():  # Ignora líneas vacías
                    parts = line.split()
                    if len(parts) == 2:
                        total_col1 += float(parts[0])
                        total_col2 += float(parts[1])

        print(f"TOTAL COLUMNA 1: {total_col1:15.2f}")
        print(f"TOTAL COLUMNA 2: {total_col2:15.2f}")

    except FileNotFoundError:
        print("Archivo DATOS.DAT no encontrado.")
    except ValueError:
        print("Error en el formato de los datos.")

if __name__ == "__main__":
    main()
