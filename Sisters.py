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
