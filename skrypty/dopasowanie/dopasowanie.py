from matplotlib import pyplot

pyplot.plotfile('dane.dat', delimiter=' ', cols=(0, 1), marker='o', color='r', ls='')

pyplot.plotfile('5/dane/delta_od_L.txt', delimiter=' ', cols=(0, 1), marker='o', color='b', ls='', markersize=3, newfig=False)
pyplot.plotfile('7/dane/delta_od_L.txt', delimiter=' ', cols=(0, 1), marker='o', color='b', ls='', markersize=3, newfig=False)
pyplot.plotfile('9/dane/delta_od_L.txt', delimiter=' ', cols=(0, 1), marker='o', color='b', ls='', markersize=3, newfig=False)
pyplot.plotfile('30/dane/delta_od_L.txt', delimiter=' ', cols=(0, 1), marker='o', color='b', ls='', markersize=3, newfig=False)

pyplot.xlim([4, 31])
pyplot.title('potencjal_chem = 0.9, masa_e = 1.0, gN0 = 0.18')

pyplot.show()