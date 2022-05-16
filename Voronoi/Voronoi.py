import matplotlib.pyplot as plt
from scipy.spatial import Voronoi, voronoi_plot_2d
points = np.array([[0, 0], [0, 1], [0, 2], [1, 0], [1, 1], [1, 2],

                   [2, 0], [2, 1], [2, 2]])
vor = Voronoi(points)
fig = voronoi_plot_2d(vor)
plt.show()
#vor.vertices
#vor.regions
#vor.ridge_points
