"""
Lancczos method on python using scipy for calculating the eigenvalues and eigenvectors.
"""
   import numpy as np
   from scipy.sparse.linalg import svds

   # Define a sparse matrix (using a NumPy array for demonstration)
   matrix = np.array([[2.0, 1.0, 1.0],
                       [1.0, 2.0, 1.0],
                       [1.0, 1.0, 2.0]])

   # Calculate a few singular values and vectors using svds
   U, s, V = svds(matrix, k=2) # k is the number of singular values/vectors to compute

   print("Singular values:", s)
   print("Left singular vectors:\n", U)
   print("Right singular vectors:\n", V)
