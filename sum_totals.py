import itertools
import os
import sys

#INDIR     = os.path.join( os.getcwd(), "Q_zeta9" )
INDIR     = os.path.join( os.getcwd(), "Q_zeta9-first_run" )
PART_SIZE = 10**4
START     = 10**0
END       = 10**8

def tabulate_class_groups( dlow, dhigh ):
	def complain(line):
		sys.stderr.write( "Got invalid class group: {oldstuff}\n".format(oldstuff=line) )
	infile_name  = os.path.join( INDIR, "{dlow}-{dhigh}.clsgps.lst".format(dlow=dlow,dhigh=dhigh) )
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

def get_disc_partitions( low, high, size ):
	n_parts = (high - low + size - 1)/size # ceil((high - low)/size)
	for i in xrange( n_parts ):
		z = min( high, low + (i+1)*size - 1 )
		yield (low + i*size, z)

partitions = get_disc_partitions( START, END, PART_SIZE )

results = {}
for partition in partitions:
	this_result = tabulate_class_groups( *partition )
	for k in this_result:
		results[k] = results.get(k,0) + this_result[k]

keys = sorted(results.keys())
for k in keys:
	print( "{class_group}\t=> {count}".format(class_group=k, count=results[k]) )

