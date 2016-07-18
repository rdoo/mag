import numpy

dic = {}

with open('Tc_od_L.txt', 'r') as file:
	for line in file:
		splitt = line.split(' ')

		if splitt[0] not in dic:
			dic[splitt[0]] = []

		dic[splitt[0]].append(float(splitt[1][:-1]))

print(dic)

with open('srednie_Tc_od_L.txt', 'w') as file:
	for item,value in dic.items():
		file.write(item + ' ' + str(numpy.mean(value)) + ' ' + str(numpy.std(value)) + '\n')