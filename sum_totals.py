import itertools
import os
import sys

if len(sys.argv) != 2:
	raise Exception("I need a list of files to process.")

FILE_LIST = None
with open(sys.argv[1]) as f:
	flines = f.readlines()
	FILE_LIST = map( lambda s: s.strip(), flines )

def tabulate_class_groups( infile_name ):
	def complain(line):
		sys.stderr.write( "Got invalid class group: {oldstuff}\n".format(oldstuff=line) )
	infile  = open( infile_name,  'r' )
	results = {}
	for line in infile:
		line = line.strip()
		# Ensure line isn't an error line
		if "ERROR" in line:
			complain(line)
			continue
		try:
			class_group_string = line.split(':')[3]
			if class_group_string == "[]":
				class_group_factors = []
			else:
				class_group_factors = map( int, str(filter(lambda c: c in "-0123456789,",class_group_string)).split(",") )
			result_key = tuple(class_group_factors)
			results[result_key] = results.get(result_key,0) + 1
		except:
			complain(line)
			continue
	infile.close()
	return results

results = {}
for infile in FILE_LIST:
	this_result = tabulate_class_groups( infile )
	for k in this_result:
		results[k] = results.get(k,0) + this_result[k]

keys = sorted(results.keys())
for k in keys:
	print( "{class_group}\t=> {count}".format(class_group=k, count=results[k]) )

