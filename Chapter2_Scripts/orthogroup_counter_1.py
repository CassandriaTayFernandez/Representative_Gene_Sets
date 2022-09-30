#Use the following python code on Orthogroups.GeneCount.tsv
#script usage:
## python orthogroup_counter.py > output.file
	x =2
	y =3
	while y < 14:
		zeroes = range(2,y)
		ones = y
		all = []
		for z in zeroes: 
			all.append("($" + str(z) + "==0)")
		all.append("($" + str(ones) + " != 0)' Orthogroups.GeneCount.tsv | wc -l")
		print(' && '.join(all))
		y += 1
