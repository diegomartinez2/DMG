#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  Sisters.py
#
#  Copyright 2026 Diego Martinez Gutierrez <diego.martinez@ehu.eus>
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
'''
A code created to solve a quaestion about who have more
sisters, boys or girls.
'''
import random
gender = ["boy", "girl"]
def generate_family(family_size):
    family_size = family_size
    family = []
    for i in range(family_size):
        family.append(random.choice(gender))
    return family
def boy_counter(family):
    boys = 0
    for sibling in family:
        if sibling == "boy":
            boys += 1
    return boys
sister_sum_for_boys = 0
boy_amount = 0
sister_sum_for_girls = 0
girl_amount = 0
for i in range(10000000):
    family = generate_family(random.randint(1, 10))
    boys = boy_counter(family)
    girls = len(family) - boys
    sister_sum_for_boys += boys*girls
    boy_amount += boys
    sister_sum_for_girls += girls*(girls-1)
    girl_amount += girls
avg_sister_for_boys = sister_sum_for_boys/boy_amount
avg_sister_for_girls = sister_sum_for_girls/girl_amount
print(avg_sister_for_girls, avg_sister_for_boys)
