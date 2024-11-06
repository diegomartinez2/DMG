#Functions to transform from lists to dictionaries and from dictionaries to lists
def list_to_dictionary(keys, values):
    """
    to_dictionary(['a', 'b'], [1, 2]) # { a: 1, b: 2 }
    """
  return dict(zip(keys, values))

def dictionary_to_list(d):
    """
    d = {'one': 1, 'three': 3, 'five': 5, 'two': 2, 'four': 4}
    dict_to_list(d)
    # [('one', 1), ('three', 3), ('five', 5), ('two', 2), ('four', 4)]
    """
  return list(d.items())
