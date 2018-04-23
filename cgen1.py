import glob
import numpy
import os
import time
import ex1_36853679 as ex1
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np

stam = ex1.fasta_iter('/home/chrx/PycharmProjects/cgen1/saccer_chrI-II.fa')
chr1 = next(stam)
chr2 = next(stam)
# ex1.test_index()
toplot = {}


results = []
for i in range(10,75):
    print(i, ',')
    for j in range(10,i):
        # print(j, ',', end='')
        seq1 = chr1[1][:i]
        seq2 = chr2[1][:j]
        s = time.time()
        ex1.align(seq1, seq2)
        e = time.time()
        tot = (e - s)
        results.append((i,j,tot))


fig = plt.figure()
ax = fig.gca(projection='3d')

# Make data.
X = np.asarray([res[0] for res in results])
Y = np.asarray([res[1] for res in results])
Z = np.asarray([res[2] for res in results])
# X = np.flip(X,0)
# Y = np.flip(Y,0)
# Z = np.flip(Z,0)
# X, Y = np.meshgrid(X, Y)
# Z = np.asarray(Z)

# x, y = X.ravel(), X.ravel()
# Plot the surface.
bottom = np.zeros_like(Z)
width = depth = 1
ax.bar3d(X, Y, bottom, width, depth, Z, shade=True)
ax.view_init(ax.elev, ax.azim-90)
plt.show()
# Customize the z axis.
# ax.set_zlim(-1.01, 1.01)
# ax.zaxis.set_major_locator(LinearLocator(10))
# ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

# Add a color bar which maps values to colors.
# fig.colorbar(surf, shrink=0.5, aspect=5)
#
# plt.show()

# for filename in glob.iglob('./reads_*.fa'):
#     print(filename)
#     try:
#         len = int(os.path.basename(filename).split('.')[0][6:])
#     except:
#         continue
#     print(len)
#     fasta_iter = ex1.fasta_iter(filename)
#     results = []
#     i = 0
#     for read in fasta_iter:
#         if i == 10:
#             break
#
#         chr1time = numpy.inf
#         chr2time = numpy.inf
#         try:
#             s = time.time()
#             ex1.align(read[1],chr1[1])
#             e = time.time()
#             chr1time = (e-s)
#         except:
#             pass
#         try:
#             s = time.time()
#             ex1.align(read[1],chr2[1])
#             e = time.time()
#             chr2time = (e-s)
#         except:
#             pass
#         results.append(min(chr1time, chr2time))
#         i+=1
#     avg = numpy.average(results)
#     toplot[len] = avg
# print(toplot)
