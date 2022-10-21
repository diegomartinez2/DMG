#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  untitled.py
#
#  Copyright 2022 Diego Martinez Gutierrez <diego.martinez@ehu.eus>
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#
#
# ---------------------------
# Importación de los módulos
# ---------------------------
import random

# -------
# Clases
# -------
class NombredeClase(object):
    """docstring for NombredeClase."""

    def __init__(self, arg):
        super(NombredeClase, self).__init__()
        self.arg = arg

# ----------
# Funciones
# ----------
def NombredeFuncion(arg):
    pass

def main(args):
    Facilitador_de_mision_list = [
    'TRABAJO DE GUANTE BLANCO PARA: Instituto de investigación de Inteligencia Artificial',
    'TRABAJO DE GUANTE BLANCO PARA: Oficina del gobierno',
    'TRABAJO DE GUANTE BLANCO PARA: Hacienda',
    'TRABAJO DE GUANTE BLANCO PARA: Industria de ciberimplantes',
    'TRABAJO DE GUANTE BLANCO PARA: Fábrica de armas',
    'TRABAJO DE GUANTE BLANCO PARA: Banco',
    'TRABAJO DE GUANTE BLANCO PARA: Planta de tratamiento de residuos',
    'TRABAJO DE GUANTE BLANCO PARA: Planta de energía',
    'TRABAJO DE GUANTE BLANCO PARA: Derechos humanos',
    'TRABAJO DE GUANTE BLANCO PARA: Fabricante de medicamentos',
    'PERSONAJE SOMBRÍO: Más cromo que humano',
    'PERSONAJE SOMBRÍO: Veterano de guerra urbana',
    'PERSONAJE SOMBRÍO: Cubierto de trapos',
    'PERSONAJE SOMBRÍO: Adicto a las drogas arruinado',
    'PERSONAJE SOMBRÍO: Arma viviente (repleto de equipo, ninja, cyberfriki)',
    'PERSONAJE SOMBRÍO: Miembro de grupo terrorista anticorporaciones',
    'PERSONAJE SOMBRÍO: Espía de megacorporación',
    'PERSONAJE SOMBRÍO: Vagabundo de las tierras baldías',
    'PERSONAJE SOMBRÍO: Mutante',
    'PERSONAJE SOMBRÍO: Esquizofrénico psicótico',
    'INTELIGENCIA ARTIFICIAL: Robot de combate',
    'INTELIGENCIA ARTIFICIAL: Androide casi humano',
    'INTELIGENCIA ARTIFICIAL: Megaordenador',
    'INTELIGENCIA ARTIFICIAL: Instalación',
    'INTELIGENCIA ARTIFICIAL: Controlado por aliens',
    'INTELIGENCIA ARTIFICIAL: De la resistencia anti-megacorpopaciones',
    'INTELIGENCIA ARTIFICIAL: Subida de la consciencia de un importante profesor fallecido',
    'INTELIGENCIA ARTIFICIAL: Autoconsciente',
    'INTELIGENCIA ARTIFICIAL: Enmascarado como lo opuesto a lo que representa',
    'INTELIGENCIA ARTIFICIAL: De origen alienigena',
    'BANDA DE ROCK: Anti corporación',
    'BANDA DE ROCK: Pro corporación',
    'BANDA DE ROCK: Consume más droga que música interpreta',
    'BANDA DE ROCK: Precursor',
    'BANDA DE ROCK: Más croma que humano',
    'BANDA DE ROCK: Subterraneo',
    'BANDA DE ROCK: Extremadamente popular',
    'BANDA DE ROCK: Más dinero y fama de la que pueden manejar',
    'BANDA DE ROCK: Banda de un solo hombre',
    'BANDA DE ROCK: Hologramas programados y computerizados y generador de música aleatorio',
    'HACKER, QUIÉN ES: 24/7 conectado y iviendo en el cyberspacio',
    'HACKER, QUIÉN ES: Miembro de un culto del ciberespacio',
    'HACKER, QUIÉN ES: Tiene rarezas, raros hábitos',
    'HACKER, QUIÉN ES: Paralizado de cuello abajo',
    'HACKER, QUIÉN ES: Abuso de droga experimental',
    'HACKER, QUIÉN ES: Famoso',
    'HACKER, QUIÉN ES: Solo un crío',
    'HACKER, QUIÉN ES: Cautivo en el ciberespacio',
    'HACKER, QUIÉN ES: Usa warez y equipo de la vieja escuela',
    'HACKER, QUIÉN ES: Operar desde el espacio',
    'StreetSecurityMilitia (SSM): Luchan contra delitos callejeros',
    'StreetSecurityMilitia (SSM): Aceptar soborno',
    'StreetSecurityMilitia (SSM): No hacen nada importante',
    'StreetSecurityMilitia (SSM): Corrupto',
    'StreetSecurityMilitia (SSM): Consta de ciudadanos preocupados',
    'StreetSecurityMilitia (SSM): Veterinarios que abusan de drogas',
    'StreetSecurityMilitia (SSM): Más cromo que humano',
    'StreetSecurityMilitia (SSM): Gatillo fácil',
    'StreetSecurityMilitia (SSM): Demasiado tímido para actuar',
    'StreetSecurityMilitia (SSM): Lucha contra megacorporación',
    'MEGACORPORACIÓN: Uno de los mayores fabricantes de medicamentos experimentales',
    'MEGACORPORACIÓN: Experimentan en humanos',
    'MEGACORPORACIÓN: Violan los derechos humanos',
    'MEGACORPORACIÓN: Precursor de I.A.',
    'MEGACORPORACIÓN: Más poder armamentístico que otra cosa',
    'MEGACORPORACIÓN: Despiadados que hacen todo lo que necesitan sin miramientos',
    'MEGACORPORACIÓN: Fabricación vehículos de lujo',
    'MEGACORPORACIÓN: Cibernética, armas, armaduras móviles...',
    'MEGACORPORACIÓN: Nadie sabe que hacen realmente',
    'MEGACORPORACIÓN: Opera desde una estación espacial',
    'DETECTIVE/OJO PRIVADO: ',
    'REPORTERO/PRESENTADOR: ',
    'MIEMBRO DE UNA BANDA: '
    ]
    Tipo_de_mision_list = [
    'PROTEGER: ',
    'INTELIGENCIA DE: ',
    'ANIQUILACIÓN: ',
    'EL HURTO: ',
    'RECUPERAR: ',
    'RESCATAR: ',
    'COMPRA-VENTA (50/50): ',
    'CAPTURAR: ',
    'OCULTAR: ',
    'PROTEGER: '
    ]
    Trabajo = random.choice(Facilitador_de_mision_list)
    Mision = random.choice(Tipo_de_mision_list)
    print (Trabajo+' quiere '+Mision)
    return 0

if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))
