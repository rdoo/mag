from math import pi

nm2au = 18.89726133921252
eV2au = 0.03674932587122423

jednostka_L = 0.286 * nm2au

with open('energia_nsp.txt', 'w') as file:
	for warstwa in range(1, 21): # max 20 warstw
		L = warstwa * jednostka_L
		file.write(str(warstwa) + ' ')
		for stan in range(1, 21): # max 20 stanow
			E = pi*pi*stan*stan/2/L/L
			file.write(str(E/eV2au) + ' ')

		file.write('\n')