from matplotlib import pyplot, rcParams, rc
from glob import glob
import numpy as numpy

rc('font', family='Arial')

rcParams['mathtext.fontset'] = 'custom'
rcParams['mathtext.rm'] = 'Bitstream Vera Sans'
rcParams['mathtext.it'] = 'Bitstream Vera Sans'
rcParams['mathtext.bf'] = 'Bitstream Vera Sans:bold'

dane_dir = 'dane/'
wykresy_dir = 'wykresy/'

#COLORS z http://tableaufriction.blogspot.ro/2012/11/finally-you-can-use-tableau-data-colors.html
niebieski = '#1F77B4'
pomaranczowy = '#FF7F0E'
zielony = '#2CA02C'
czerwony = '#D62728'
fioletowy = '#9467BD'


def delta_od_L():
	pyplot.plotfile(dane_dir + 'delta_od_L.txt', delimiter=' ', names=['a', 'b'], cols=(0, 1), marker='.', color=niebieski, ls='')
	pyplot.plot((1, 5), (0.4, 0.4), 'k-')
	pyplot.text(1, 0.5, 'wartosc bulk')
	pyplot.xlabel('L [nm]')
	pyplot.ylabel(u'\u0394 [eV]')

	pyplot.savefig(wykresy_dir + 'delta_od_L.png')


def Tc_od_L():
	pyplot.plotfile(dane_dir + 'Tc_od_L.txt', delimiter=' ', names=['a', 'b'], cols=(0, 1), marker='.', color=niebieski, ls='')
	pyplot.plot((1, 5), (4, 4), 'k-')
	pyplot.text(1, 4, 'wartosc bulk')
	pyplot.xlabel('L [nm]')
	pyplot.ylabel(r'$T_c$ [K]')

	pyplot.savefig(wykresy_dir + 'Tc_od_L.png')

def potencjal_od_L():
	pyplot.plotfile(dane_dir + 'potencjal_od_L.txt', delimiter=' ', names=['a', 'b'], cols=(0, 1), marker='.', color=niebieski, ls='')
	pyplot.plot((1, 5), (0.9, 0.9), 'k-')
	pyplot.text(1, 0.95, 'wartosc bulk')
	pyplot.xlabel('L [nm]')
	pyplot.ylabel(u'\u03bc [eV]')

	pyplot.savefig(wykresy_dir + 'potencjal_od_L.png')

def delta_od_z():
	pliki = glob(dane_dir + 'delta_od_z*.txt')

	for nazwa in pliki:
			pyplot.plotfile(nazwa, delimiter=' ', names=['a', 'b'], cols=(0, 1), marker='.', color=niebieski, ls='')
			pyplot.xlabel('z [nm]')
			pyplot.ylabel(u'\u03bc [eV]') # czy mev?

			pyplot.savefig(nazwa.replace('dane', 'wykresy').replace('txt', 'png'))

def delta_od_T():
	pliki = glob(dane_dir + 'delta_od_T*.txt')

	for nazwa in pliki:
			pyplot.plotfile(nazwa, delimiter=' ', names=['a', 'b'], cols=(0, 1), marker='.', color=niebieski, ls='')
			pyplot.xlabel('T [K]')
			pyplot.ylabel(u'\u03bc [eV]') # czy mev?

			pyplot.savefig(nazwa.replace('dane', 'wykresy').replace('txt', 'png'))

def dyspersja():
	pliki = glob(dane_dir + 'dyspersja*.txt')

	for nazwa in pliki:
		k, e1, e2, e3, e4, e5, e6, e7, e8 = numpy.loadtxt(nazwa, unpack=True)

		pyplot.clf()

		pyplot.plot(k, e1, color=niebieski)
		pyplot.plot(k, e2, color=pomaranczowy)
		pyplot.plot(k, e3, color=zielony)
		pyplot.plot(k, e4, color=fioletowy)

		pyplot.xlabel('k [1/nm]')
		pyplot.ylabel('E [meV]')

		new_axes = pyplot.gcf().add_axes([0.2, 0.6, 0.25, 0.25])

		new_axes.plot(k, e5, color=niebieski)
		new_axes.plot(k, e6, color=pomaranczowy)
		new_axes.plot(k, e7, color=zielony)
		new_axes.plot(k, e8, color=fioletowy)

		pyplot.xlabel('k [1/nm]')
		pyplot.ylabel(r'$E_kin$ [meV]')

		pyplot.savefig(nazwa.replace('dane\\', 'wykresy\\t').replace('txt', 'png'))



		pyplot.clf()

		pyplot.plot(k, e1, color=niebieski)
		pyplot.plot(k, e2, color=pomaranczowy)
		pyplot.plot(k, e3, color=zielony)
		pyplot.plot(k, e4, color=fioletowy)

		pyplot.xlabel('k [1/nm]')
		pyplot.ylabel('E [meV]')
		pyplot.ylim(0, 10)

		new_axes = pyplot.gcf().add_axes([0.2, 0.6, 0.25, 0.25])

		new_axes.plot(k, e5, color=niebieski)
		new_axes.plot(k, e6, color=pomaranczowy)
		new_axes.plot(k, e7, color=zielony)
		new_axes.plot(k, e8, color=fioletowy)

		pyplot.xlabel('k [1/nm]')
		pyplot.ylabel(r'$E_{kin}$ [meV]')
		pyplot.ylim(-50, 50)

		pyplot.savefig(nazwa.replace('dane', 'wykresy').replace('txt', 'png'))


delta_od_L()
Tc_od_L()
potencjal_od_L()
delta_od_z()
delta_od_T()
dyspersja()