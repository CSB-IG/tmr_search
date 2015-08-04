# This script helps to collapse the archives of microarray expression
# It is required that the first and second colunma are the Gene Symbol and the statistical b consecutively ("loads" in our case)
# It runs like any Python script:
# $ python microarray_colapser.py path/to/the/exp_matrix_file.txt
# output in this example would exp_matrix_file_colapsed.txt
# The output file still has the original first two columns


import sys



exp_file = open(sys.argv[1], 'rb')

new_lines = {}
dict1 = {}
names = []

for line in exp_file :
	p = line.split("\t")
	if p[0] in dict1:
		if p[1] > dict1[p[0]]:
			dict1 [p[0]] = p[1]
			new_lines [p[0]] = line.strip()
	else:
		dict1 [p[0]] = p[1]
		new_lines [p[0]] = line.strip()
		names.append(p[0])
		
exp_file.close()

output = open(sys.argv[1].replace(".txt", "", 1)+"_colapsed.txt", "w")
for n in names:
	output.write(new_lines[n]+"\n")
output.close()	