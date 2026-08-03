def texto_a_braille(texto):
    # Diccionario de traducción simplificado (Braille español / internacional)
    equivalencias_braille = {
        'a': '⠁', 'b': '⠃', 'c': '⠉', 'd': '⠙', 'e': '⠑', 'f': '⠋', 'g': '⠛',
        'h': '⠃', 'i': '⠊', 'j': '⠚', 'k': '⠁', 'l': '⠇', 'm': '⠍', 'n': '⠝',
        'ñ': '⠌', 'o': '⠕', 'p': '⠏', 'q': '⠟', 'r': '⠌', 's': '⠎', 't': '⠞',
        'u': '⠥', 'v': '⠧', 'w': '⠺', 'x': '⠭', 'y': '⠽', 'z': '⠵',
        'á': '⠷', 'é': '⠮', 'í': '⠌', 'ó': '⠹', 'ú': '⠾',
        ' ': ' ', '.': '⠲', ',': '⠂', '!': '⠔', '?': '⠢'
    }
    
    resultado = []
    for caracter in texto.lower():
        # Busca el símbolo braille, si no existe lo deja igual
        resultado.append(equivalencias_braille.get(caracter, caracter))
        
    return "".join(resultado)

# Prueba del traductor
mensaje = "hola mundo"
print(f"Texto original: {mensaje}")
print(f"En Braille:    {texto_a_braille(mensaje)}")
