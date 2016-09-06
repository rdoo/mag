from matplotlib import pyplot, rcParams, rc

rc('font', family='Arial')

rcParams['mathtext.fontset'] = 'custom'
rcParams['mathtext.rm'] = 'Bitstream Vera Sans'
rcParams['mathtext.it'] = 'Bitstream Vera Sans'
rcParams['mathtext.bf'] = 'Bitstream Vera Sans:bold'

niebieski = '#1F77B4'

pyplot.plotfile('1/dane/delta_od_L.txt', delimiter=' ', names=['a', 'b'], cols=(0, 1), marker='.', color=niebieski, ls='')
pyplot.plotfile('2/dane/delta_od_L.txt', delimiter=' ', names=['a', 'b'], cols=(0, 1), marker='.', color=niebieski, ls='', newfig=False)
pyplot.plotfile('3/dane/delta_od_L.txt', delimiter=' ', names=['a', 'b'], cols=(0, 1), marker='.', color=niebieski, ls='', newfig=False)
pyplot.plotfile('4/dane/delta_od_L.txt', delimiter=' ', names=['a', 'b'], cols=(0, 1), marker='.', color=niebieski, ls='', newfig=False)

pyplot.xlabel('L [nm]')
pyplot.ylabel(u'\u0394 [eV]')

pyplot.savefig('delta_od_L.png')


pyplot.clf()

pyplot.plotfile('1/dane/Tc_od_L.txt', delimiter=' ', names=['a', 'b'], cols=(0, 1), marker='.', color=niebieski, ls='')
pyplot.plotfile('2/dane/Tc_od_L.txt', delimiter=' ', names=['a', 'b'], cols=(0, 1), marker='.', color=niebieski, ls='', newfig=False)
pyplot.plotfile('3/dane/Tc_od_L.txt', delimiter=' ', names=['a', 'b'], cols=(0, 1), marker='.', color=niebieski, ls='', newfig=False)
pyplot.plotfile('4/dane/Tc_od_L.txt', delimiter=' ', names=['a', 'b'], cols=(0, 1), marker='.', color=niebieski, ls='', newfig=False)

pyplot.xlabel('L [nm]')
pyplot.ylabel(r'$T_c$ [K]')

pyplot.savefig('Tc_od_L.png')