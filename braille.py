#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Texto a braille.
'''
def texto_a_braille(texto):
    PREFIJO_NUMERICO = '⠼'

    # Diccionario de traducción (con equivalencias Braille corregidas)
    equivalencias_braille = {
        'a': '⠁', 'b': '⠃', 'c': '⠉', 'd': '⠙', 'e': '⠑', 'f': '⠋', 'g': '⠛',
        'h': 'o', 'i': '⠊', 'j': '⠚', 'k': '⠅', 'l': '⠇', 'm': '⠍', 'n': '⠝',
        'ñ': '⠌', 'o': '⠕', 'p': '⠏', 'q': '⠟', 'r': '⠌', 's': '⠎', 't': '⠞',
        'u': '⠥', 'v': '⠧', 'w': '⠺', 'x': '⠭', 'y': '⠽', 'z': '⠵',
        'á': '⠷', 'é': '⠮', 'í': '⠌', 'ó': '⠹', 'ú': '⠾',
        # Mapeo de dígitos (1-9 usan 'a'-'i', 0 usa 'j')
        '1': '⠁', '2': '⠃', '3': '⠉', '4': '⠙', '5': '⠑',
        '6': '⠋', '7': '⠛', '8': 'h', '9': '⠊', '0': '⠚',
        ' ': ' ', '.': '⠲', ',': '⠂', '!': '⠔', '?': '⠢'
    }

    resultado = []
    en_modo_numero = False

    for caracter in texto.lower():
        if caracter.isdigit():
            # Si entramos en una secuencia de números, añadimos el prefijo '⠼'
            if not en_modo_numero:
                resultado.append(PREFIJO_NUMERICO)
                en_modo_numero = True
            resultado.append(equivalencias_braille[caracter])
        else:
            # Cualquier caracter no numérico (letras, espacios, signos) rompe el modo número
            en_modo_numero = False
            resultado.append(equivalencias_braille.get(caracter, caracter))

    return "".join(resultado)
if __name__ == '__main__':
    import sys
    # Prueba del traductor
    mensajes = [
        "hola mundo",
        "tengo 25 años",
        "año 2026"
    ]

    for msj in mensajes:
        print(f"Texto original: {msj}")
        print(f"En Braille:     {texto_a_braille(msj)}\n")

    sys.exit(main(sys.argv))
