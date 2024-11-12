# This function merge lists in python.
# from [a,b],[1,2] to [[a,1],[b,2]]
def merge(*args, fill_value = None):
    """
        Merge lists in python, examples:
    merge(['a', 'b'], [1, 2], [True, False]) # [['a', 1, True], ['b', 2, False]]
    merge(['a'], [1, 2], [True, False]) # [['a', 1, True], [None, 2, False]]
    merge(['a'], [1, 2], [True, False], fill_value = '_')
    # [['a', 1, True], ['_', 2, False]]
    """
  max_length = max([len(lst) for lst in args])
  result = []
  for i in range(max_length):
    result.append([
      args[k][i] if i < len(args[k]) else fill_value for k in range(len(args))
    ])
  return result


# Simple sorting function for lists

def sort_by_indexes(lst, indexes, reverse=False):
    """
        Sort lists by index, examples:
a = ['eggs', 'bread', 'oranges', 'jam', 'apples', 'milk']
b = [3, 2, 6, 4, 1, 5]
sort_by_indexes(a, b) # ['apples', 'bread', 'eggs', 'jam', 'milk', 'oranges']
sort_by_indexes(a, b, True)
# ['oranges', 'milk', 'jam', 'eggs', 'bread', 'apples']
    """
  return [val for (_, val) in sorted(zip(indexes, lst), key=lambda x: \
          x[0], reverse=reverse)]
