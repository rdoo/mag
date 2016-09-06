from math import pi
from matplotlib import pyplot, rcParams, rc
from glob import glob
import numpy as numpy

rc('font', family='Arial')

rcParams['mathtext.fontset'] = 'custom'
rcParams['mathtext.rm'] = 'Bitstream Vera Sans'
rcParams['mathtext.it'] = 'Bitstream Vera Sans'
rcParams['mathtext.bf'] = 'Bitstream Vera Sans:bold'

niebieski = '#1F77B4'

nm2au = 18.89726133921252
eV2au = 0.03674932587122423

jednostka_L = 0.286 * nm2au

l1, p1 = numpy.loadtxt('1/dane/potencjal_od_L.txt', unpack=True)
l1, p2 = numpy.loadtxt('2/dane/potencjal_od_L.txt', unpack=True)
l1, p3 = numpy.loadtxt('3/dane/potencjal_od_L.txt', unpack=True)
l1, p4 = numpy.loadtxt('4/dane/potencjal_od_L.txt', unpack=True)

tablica_mi = p1.tolist() + p2.tolist() + p3.tolist() + p4.tolist()

print(tablica_mi)

# tablica_mi = [1.99951171875000000000, 1.68017578125000000000, 1.30712890625000000000, 1.24951171875000022204, 1.12353515625000000000, 1.11572265625000000000, 1.05712890624999977796, 1.04345703125000022204, 1.02294921875000000000, 1.00830078125000000000, 1.00244140624999977796, 0.98876953125000022204, 0.98779296875000011102, 0.97509765624999988898, 0.97314453125000000000, 0.96630859375000022204, 0.96240234374999988898, 0.95947265625000011102, 0.95556640624999977796, 0.95361328125000000000]

with open('energia_nsp.txt', 'w') as file:
	for warstwa in range(1, 21): # max 20 warstw
		L = warstwa * jednostka_L
		file.write(str(warstwa) + ' ')
		for stan in range(1, 21): # max 20 stanow
			E = pi*pi*stan*stan/2/L/L - tablica_mi[warstwa-1] * eV2au
			file.write(str(E/eV2au) + ' ')

		file.write('\n')

Ls = []
Es = []

for warstwa in range(1, 21):
	L = warstwa * jednostka_L
	for stan in range(1, 21): # max 20 stanow
		E = pi*pi*stan*stan/2/L/L - tablica_mi[warstwa-1] * eV2au
		Ls.append(warstwa)
		Es.append(E/eV2au)

pyplot.plot(Ls, Es, marker='.', ls='', color=niebieski)
pyplot.plot((0, 20), (0, 0), 'k-')
pyplot.ylim(-1, 1)

pyplot.savefig('energia_nsp.png')