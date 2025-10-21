import numpy as np
from scipy.linalg import eig
import os

# =======================================================
# FUNCIONES AUXILIARES (Traducción de Funciones FORTRAN)
# =======================================================

def BIFACT(mm: int) -> int:
    """
    Calcula el doble factorial (!!).
    (mm)!! = mm * (mm-2) * (mm-4) * ...
    """
    if mm < -1:
        # El doble factorial de números impares negativos es complejo y no se usa aquí.
        # Basado en la implementación de FORTRAN, solo se manejan >= 0 o -1.
        # Si mm es par y negativo, la función FORTRAN devuelve 1, pero debe ser indefinido.
        # Para evitar problemas, solo trabajamos con positivos, como en el contexto de F().
        return 0

    if mm <= 0:
        return 1 # (0)!! = 1 y (-1)!! = 1

    result = 1
    for i in range(mm, 0, -2):
        result *= i
    return result

def F(a: int, b: int, TIPO: int, alpha: float = 0.0) -> float:
    """
    Función de integral F(a, b, TIPO).
    TIPO=1 (cuadrado), TIPO=2 (cilíndrico), TIPO=3 (hexagonal).
    """
    # La función F en FORTRAN usa la propiedad de que la integral es cero si a o b es impar.
    # El código se simplifica usando (1 + (-1)**k) que es 2 si k es par y 0 si k es impar.

    PI = np.pi

    if TIPO == 1:
        # Forma cuadrada: F=(1.d0+(-1.d0)**a)/(2.d0*(a+1))*(1.d0+(-1.d0)**b)/(2.d0*(b+1))
        term_a = (1.0 + (-1.0)**a) / (2.0 * (a + 1))
        term_b = (1.0 + (-1.0)**b) / (2.0 * (b + 1))
        return term_a * term_b

    elif TIPO == 2:
        # Forma cilíndrica (Circular):
        # F=(1.d0-alpha**(a+b+2))*(4.d0*PI*BIFACT(a-1)*BIFACT(b-1)/BIFACT(a+b+2))*((1.d0+(-1.d0)**a)*(1.d0+(-1.d0)**b))/4.d0
        if (a % 2 != 0) or (b % 2 != 0):
            return 0.0

        term_bifact = (4.0 * PI * BIFACT(a - 1) * BIFACT(b - 1)) / BIFACT(a + b + 2)
        term_alpha = (1.0 - alpha**(a + b + 2))
        # La parte final ((1+(-1)**a)*(1+(-1)**b))/4.d0 es 1 si a y b son pares, 0 si no.
        # Ya cubierto por el chequeo de paridad.
        return term_alpha * term_bifact

    elif TIPO == 3:
        # Forma hexagonal:
        # F=((1.d0-alpha**2)/((2**(a-1))*(b+1)*(a+b+2)))*((1.d0+(-1.d0)**a)*(1.d0+(-1.d0)**b))/4.d0
        if (a % 2 != 0) or (b % 2 != 0):
            return 0.0

        term_geom = (1.0 - alpha**2) / ((2**(a - 1)) * (b + 1) * (a + b + 2))
        return term_geom

    else:
        raise ValueError("Error: TIPO de sección transversal no reconocido.")


# =======================================================
# PROGRAMA PRINCIPAL
# =======================================================

def infwire_calculator():
    """Ejecuta el programa infwire para calcular modos normales acústicos."""

    # --- 1. CONFIGURACIÓN Y CONSTANTES ---
    N = 12
    # NDIM = (N * (N + 1)) / 2  * 3 (por el índice espacial)
    # La lógica de FORTRAN es (número de m,n) * 3 direcciones.
    # El bucle de FORTRAN cuenta N_TERMS = (N+1)(N+2)/2, pero filtra con m+n <= N.
    # El valor máximo de INDICE_GENERAL es (N+1)(N+2)/2 * 3.
    # Para N=12, (13*14)/2 = 91. NDIM = 91 * 3 = 273.
    # La declaración 'parameter (NDIM=(N*(N+1))/2)' en FORTRAN parece ser un error de copia/pega,
    # ya que después el programa cuenta INDICE_GENERAL hasta 273 (para N=12) y falla.
    # Corregiremos esto usando el INDICE_GENERAL calculado en tiempo de ejecución.

    # Definimos NDIM con un valor seguro. Lo ajustaremos después de contar.
    NDIM_MAX = int(3 * (N + 1) * (N + 2) / 2)

    # UI = Numero 'i' (complex*16)
    UI = 1j

    # Matriz R(i, j) para relacionar C(i,j,k,l) con C(a,b)
    # R(indice_i, 1) -> R(i,1) = 1 (xx), 6 (yy), 5 (zz) -> ERROR: La matriz R de FORTRAN es 3x3
    # R: [1, 6, 5], [6, 2, 4], [5, 4, 3]
    # R[i, j] mapea el índice de desplazamiento (1=x, 2=y, 3=z) a la notación de Voigt (1..6)
    R = np.array([
        [1, 6, 5],  # i=1 (x): C(x,x,...) -> C(1,...); C(x,y,...) -> C(6,...); C(x,z,...) -> C(5,...)
        [6, 2, 4],  # i=2 (y): C(y,x,...) -> C(6,...); C(y,y,...) -> C(2,...); C(y,z,...) -> C(4,...)
        [5, 4, 3]   # i=3 (z): C(z,x,...) -> C(5,...); C(z,y,...) -> C(4,...); C(z,z,...) -> C(3,...)
    ], dtype=int) - 1 # Restamos 1 para indexación de Python (0-based)

    # --- 2. LECTURA DE DATOS ---
    try:
        with open('hilo_inf.dat', 'r') as f:
            Rho, dim_A, dim_B, TIPO = map(float, f.readline().split())
            Rho, TIPO = float(Rho), int(TIPO)
    except FileNotFoundError:
        print("ADVERTENCIA: Archivo 'hilo_inf.dat' no encontrado. Usando valores por defecto.")
        Rho, dim_A, dim_B, TIPO = 1000.0, 0.5, 0.5, 1

    try:
        C_flat = np.loadtxt('c.dat')
        C = C_flat.reshape(6, 6)
    except FileNotFoundError:
        print("ADVERTENCIA: Archivo 'c.dat' no encontrado. Usando matriz C isotrópica por defecto.")
        C11, C12, C44 = 280e9, 120e9, 80e9 # MPa -> Pa para unidades consistentes
        C = np.zeros((6, 6))
        C[0, 0] = C[1, 1] = C[2, 2] = C11
        C[0, 1] = C[0, 2] = C[1, 0] = C[1, 2] = C[2, 0] = C[2, 1] = C12
        C[3, 3] = C[4, 4] = C[5, 5] = C44

    # --- 3. GENERACIÓN DE ÍNDICES DE MODOS ---

    # Variables de almacenamiento de índices (dimensiones max NDIM_MAX)
    indice_espacial = np.zeros(NDIM_MAX, dtype=int)
    etiqueta_m = np.zeros(NDIM_MAX, dtype=int)
    etiqueta_n = np.zeros(NDIM_MAX, dtype=int)

    INDICE_GENERAL = 0

    for i in range(1, 4):  # indice espacial: 1 (x), 2 (y), 3 (z)
        for m in range(0, N + 1):
            for ne in range(0, N + 1):
                suma = m + ne
                if suma <= N:
                    INDICE_GENERAL += 1
                    indice_espacial[INDICE_GENERAL - 1] = i
                    etiqueta_m[INDICE_GENERAL - 1] = m
                    etiqueta_n[INDICE_GENERAL - 1] = ne

    NDIM = INDICE_GENERAL
    print(f"Número de términos de la expansión (NDIM): {NDIM}")

    # Recortar arrays a la dimensión real
    indice_espacial = indice_espacial[:NDIM]
    etiqueta_m = etiqueta_m[:NDIM]
    etiqueta_n = etiqueta_n[:NDIM]

    # --- 4. CÁLCULO DE MATRICES Gamma (Compleja) y E (Real) ---

    # Bucle 'do j=0,0' -> solo q=0
    q = 0.0

    # Inicialización de matrices
    Gamma = np.zeros((NDIM, NDIM), dtype=complex)
    E = np.zeros((NDIM, NDIM), dtype=float)

    print("Calculando matrices Gamma y E...")

    for ii in range(NDIM):
        for jj in range(ii, NDIM): # Solo calculamos la mitad superior (matrices simétricas)

            # Índices 0-based de Python
            indice_i = indice_espacial[ii] - 1 # 0, 1, 2
            indice_j = indice_espacial[jj] - 1 # 0, 1, 2

            m_i = etiqueta_m[ii]
            n_i = etiqueta_n[ii]
            m_j = etiqueta_m[jj]
            n_j = etiqueta_n[jj]

            # La fórmula de Gamma tiene 9 términos (producto de 3x3 tensiones)
            gamma_val = 0.0 + 0.0j

            # 1. C(R(i,1), R(j,1)) * m_i * m_j * F(m_i + m_j - 2, n_i + n_j)
            gamma_val += C[R[indice_i, 0], R[indice_j, 0]] * m_i * m_j * \
                         F(m_i + m_j - 2, n_i + n_j, TIPO)

            # 2. C(R(i,2), R(j,2)) * n_i * n_j * F(m_i + m_j, n_i + n_j - 2)
            gamma_val += C[R[indice_i, 1], R[indice_j, 1]] * n_i * n_j * \
                         F(m_i + m_j, n_i + n_j - 2, TIPO)

            # 3. C(R(i,3), R(j,3)) * Q * Q * F(m_i + m_j, n_i + n_j)
            # Como q=0, este término es 0.
            # gamma_val += C[R[indice_i, 2], R[indice_j, 2]] * q * q * \
            #              F(m_i + m_j, n_i + n_j, TIPO)

            # 4. C(R(i,1), R(j,2)) * m_i * n_j * F(m_i + m_j - 1, n_i + n_j - 1)
            gamma_val += C[R[indice_i, 0], R[indice_j, 1]] * m_i * n_j * \
                         F(m_i + m_j - 1, n_i + n_j - 1, TIPO)

            # 5. C(R(i,1), R(j,3)) * m_i * UI * Q * F(...)
            # Como q=0, este término es 0.
            # gamma_val += C[R[indice_i, 0], R[indice_j, 2]] * m_i * UI * q * \
            #              F(m_i + m_j - 1, n_i + n_j, TIPO)

            # 6. C(R(i,2), R(j,1)) * n_i * m_j * F(m_i + m_j - 1, n_i + n_j - 1)
            gamma_val += C[R[indice_i, 1], R[indice_j, 0]] * n_i * m_j * \
                         F(m_i + m_j - 1, n_i + n_j - 1, TIPO)

            # 7. C(R(i,2), R(j,3)) * n_i * UI * Q * F(...)
            # Como q=0, este término es 0.
            # gamma_val += C[R[indice_i, 1], R[indice_j, 2]] * n_i * UI * q * \
            #              F(m_i + m_j - 1, n_i + n_j, TIPO)

            # 8. C(R(i,3), R(j,1)) * (-1) * UI * Q * m_j * F(...)
            # Como q=0, este término es 0.
            # gamma_val += C[R[indice_i, 2], R[indice_j, 0]] * (-1) * UI * q * m_j * \
            #              F(m_i + m_j - 1, n_i + n_j, TIPO)

            # 9. C(R(i,3), R(j,2)) * (-1) * UI * Q * n_j * F(...)
            # Como q=0, este término es 0.
            # gamma_val += C[R[indice_i, 2], R[indice_j, 1]] * (-1) * UI * q * n_j * \
            #              F(m_i + m_j, n_i + n_j - 1, TIPO)

            # Asignar a la matriz Gamma (Gamma(ii,jj) = Gamma(jj,ii) porque es simétrica en C, m, n, F)
            Gamma[ii, jj] = gamma_val
            Gamma[jj, ii] = gamma_val

            # Cálculo de la matriz E (diagonal por bloque)
            E_val = 0.0
            if indice_i == indice_j: # Solo cuando los índices espaciales son iguales (x-x, y-y, z-z)
                E_val = Rho * F(etiqueta_m[ii] + etiqueta_m[jj],
                                etiqueta_n[ii] + etiqueta_n[jj], TIPO)

            E[ii, jj] = E_val
            E[jj, ii] = E_val


    # --- 5. SOLUCIÓN DEL PROBLEMA DE VALORES PROPIOS ---

    # El problema es: Gamma * v = omega^2 * E * v
    # SciPy: eig(A, B) resuelve A * v = lambda * B * v, donde lambda = omega^2.

    # Notar que en FORTRAN se resuelve E*v = lambda*Gamma*v, donde lambda = beta/alpha
    # Para ser estrictos con la notación, resolveremos Gamma * v = lambda * E * v
    print("Resolviendo el problema generalizado de valores propios...")

    # Los autovalores (w_sq) son lambda = omega^2.
    # Los autovectores (vects) son v.
    w_sq, vects = eig(Gamma, E)

    # --- 6. ESCRITURA DE RESULTADOS ---

    # El FORTRAN guarda sqrt(omega) / (2*PI)
    # Omega = sqrt(lambda)
    # Frecuencia (Hz) = Omega / (2*PI)

    # Filtramos valores espurios (parte imaginaria grande)
    w_sq_real = np.real(w_sq)

    # Solo consideramos autovalores no negativos, ya que omega^2 debe ser real y positivo.
    # Los valores negativos o altamente imaginarios son numéricos o modos evanescentes.
    eigenfrequencies = []

    for i, w2 in enumerate(w_sq):
        if w2.real >= 0 and abs(w2.imag) < 1e-9 * w2.real: # Ignorar la pequeña parte imaginaria
            omega_val = np.sqrt(w2.real)
            frecuencia = omega_val / (2.0 * np.pi)
            eigenfrequencies.append((i + 1, frecuencia))
        # elif w2.real < 0 and abs(w2.imag) < 1e-9 * abs(w2.real):
        #     # Modos evanescentes puros
        #     frecuencia_imag = np.sqrt(-w2.real) / (2.0 * np.pi)
        #     eigenfrequencies.append((i + 1, 1j * frecuencia_imag))
        else:
            # Modos complejos o negativos (generalmente atenuantes o numéricos)
            pass

    # Se imprimen los valores propios en un archivo
    output_file = 'salida_frecuencias_python.dat'
    print(f"Escribiendo {len(eigenfrequencies)} frecuencias reales a '{output_file}'...")

    with open(output_file, 'w') as f:
        f.write("# ID | Frecuencia (Hz)\n")
        for i, freq in eigenfrequencies:
            f.write(f"{i}\t{freq:.12e}\n")

    print("\nCálculo finalizado. Revise 'salida_frecuencias_python.dat'.")

if __name__ == "__main__":
    infwire_calculator()
