from matplotlib import pyplot

pyplot.plotfile('delta_od_L1.txt', delimiter=' ', cols=(0, 1), marker='o', color='b', ls='')
pyplot.plotfile('delta_od_L2.txt', delimiter=' ', cols=(0, 1), marker='o', color='b', ls='', newfig=False)
pyplot.plotfile('delta_od_L3.txt', delimiter=' ', cols=(0, 1), marker='o', color='b', ls='', newfig=False)
pyplot.plotfile('delta_od_L4.txt', delimiter=' ', cols=(0, 1), marker='o', color='b', ls='', newfig=False)


pyplot.title('potencjal_chem = 0.9, masa_e = 1.0, gN0 = 0.18')

pyplot.show()