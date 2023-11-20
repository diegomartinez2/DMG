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
