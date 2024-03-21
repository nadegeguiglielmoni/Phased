#/usr/bin/python3

import Bio
from Bio import SeqIO
import argparse

def parse_args():
	""" 
	Gets the arguments from the command line.
	"""

	parser = argparse.ArgumentParser()
	parser.add_argument('-m', '--mapping', required=True,
				help="""Mapped long reads.""")
	parser.add_argument('-d', '--dups', required=True,
				help="""Haplotig matches.""")
	return parser.parse_args()


def main():

	args = parse_args()
	mapping = args.mapping
	dups = args.dups

	infile = open(dups, "r")

	haplotigs = {}

	for line in infile.readlines():
		line = line.strip().split()
		if "HAPLOTIG" in line[3]:
			if line[0] not in haplotigs.keys():
				haplotigs[line[0]] = {}
			haplotigs[line[0]][line[4]] = []

	#print(haplotigs)

	infile = open(mapping, "r")

	count_reads = {}
	reads_coord = {}

	for line in infile.readlines():
		line = line.strip().split()
		if line[0] not in reads_coord.keys():
			coord = (str(line[5]), int(line[2]), int(line[3]))
			reads_coord[str(line[0])] = [coord]
		else :
			coord = (str(line[5]), int(line[2]), int(line[3]))
			reads_coord[str(line[0])].append(coord)

	#print(reads_coord)

	print(str(len(reads_coord.keys())) + " reads in total")

	for read in reads_coord.keys():
		if read not in count_reads.keys():
			count_reads[read] = 0
		for i in range(0, len(reads_coord[read])):
			for j in range(1, len(reads_coord[read])):
				#print(read, reads_coord[read])
				if reads_coord[read][i][0] is not reads_coord[read][j][0]:
					if reads_coord[read][i][1] < reads_coord[read][j][1] and reads_coord[read][i][2] < reads_coord[read][j][2]:
						#print(reads_coord[read][i][0],reads_coord[read][j][0])
						if reads_coord[read][i][0] in haplotigs.keys() and reads_coord[read][j][0] in haplotigs[reads_coord[read][i][0]].keys() and reads_coord[read][i][2] <= reads_coord[read][j][1]:
							if read not in haplotigs[reads_coord[read][i][0]][reads_coord[read][j][0]]:
								haplotigs[reads_coord[read][i][0]][reads_coord[read][j][0]].append(read)
						elif reads_coord[read][j][0] in haplotigs.keys() and j in haplotigs[reads_coord[read][j][0]].keys():
							if read not in haplotigs[reads_coord[read][i][0]][reads_coord[read][j][0]]:
								haplotigs[reads_coord[read][j][0]][reads_coord[read][i][0]].append(read)
					elif reads_coord[read][i][1] > reads_coord[read][j][1] and reads_coord[read][i][2] > reads_coord[read][j][2]:
						if reads_coord[read][i][0] in haplotigs.keys() and reads_coord[read][j][0] in haplotigs[reads_coord[read][i][0]].keys() and reads_coord[read][i][1] >= reads_coord[read][j][2]:
							if read not in haplotigs[reads_coord[read][i][0]][reads_coord[read][j][0]]:
								haplotigs[reads_coord[read][i][0]][reads_coord[read][j][0]].append(read)
						elif reads_coord[read][j][0] in haplotigs.keys() and j in haplotigs[reads_coord[read][j][0]].keys():
							if read not in haplotigs[reads_coord[read][i][0]][reads_coord[read][j][0]]:
								haplotigs[reads_coord[read][j][0]][reads_coord[read][i][0]].append(read)

	#print(haplotigs)

	count = 0

	for i in haplotigs.keys():
		for j in haplotigs[i].keys(): 
			#print(i, j, len(haplotigs[i][j]))
			count = count + len(haplotigs[i][j])

	print(str(count) + " switched reads")

	print(count/len(reads_coord.keys())*100)

if __name__ == '__main__':
	main()

