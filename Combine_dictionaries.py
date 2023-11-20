from collections import defaultdict

def combine_values(*dicts):
    """
    Combine dictionaries into one, example:
    d1 = {'a': 1, 'b': 'foo', 'c': 400}
    d2 = {'a': 3, 'b': 200, 'd': 400}

    combine_values(d1, d2) # {'a': [1, 3], 'b': ['foo', 200], 'c': [400], 'd': [400]}
    """
  res = defaultdict(list)
  for d in dicts:
    for key in d:
      res[key].append(d[key])
  return dict(res)

from collections import defaultdict

def collect_dictionary(obj):
    """
    Inverts a dictionary with non-unique hashable values. Example:
ages = {
  'Peter': 10,
  'Isabel': 10,
  'Anna': 9,
}
collect_dictionary(ages) # { 10: ['Peter', 'Isabel'], 9: ['Anna'] }    
    """
  inv_obj = defaultdict(list)
  for key, value in obj.items():
    inv_obj[value].append(key)
  return dict(inv_obj)
