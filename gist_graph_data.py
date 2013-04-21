import itertools
import os
import sys

INDIR     = os.path.join( os.getcwd(), "Q_zeta9" )
OUTDIR    = os.path.join( os.getcwd(), "Q_zeta9-analysis" )
OUTFILE_T = os.path.join( OUTDIR, "graph_p={prime}_{class_gp}.dat" )
PART_SIZE = 10**4
START     = 10**0
END       = 10**8
PRIMES    = [3, 5, 7]
CLASS_GPS = \
	{ 3: [(), (3,), (3,3), (9,)]
	, 5: [(), (5,), (5,5), (25,)]
	, 7: [(), (7,), (7,7), (49,)]
	}

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

def reduce_to_p_part( results, p ):
	new_results = {}
	for decomposition in results:
		ppart = filter( lambda n: n % p == 0, decomposition )
		ppart = tuple(ppart)
		new_results[ppart] = new_results.get(ppart,0) + results[decomposition]
	return new_results

partitions = get_disc_partitions( START, END, PART_SIZE )

files = {}
for p in PRIMES:
	files[p] = {}
	for class_gp in CLASS_GPS[p]:
		if class_gp == ():
			class_gp_name = "trivial"
		else:
			class_gp_name = str(class_gp[0])
			for e in class_gp[1:]:
				class_gp_name += "," + str(e)
		files[p][class_gp] = open( OUTFILE_T.format(prime=p, class_gp=class_gp_name), 'w' )

results = {}
for partition in partitions:
	this_result = tabulate_class_groups( *partition )
	for k in this_result:
		results[k] = results.get(k,0) + this_result[k]
	for p in PRIMES:
		ppart = reduce_to_p_part(results, p)
		total = 0
		for k in ppart:
			total += ppart[k]
		for class_gp in CLASS_GPS[p]:
			this_stat = ppart.get(class_gp,0)
			percent = 100.0*float(this_stat)/float(total)
			outline = "{disc_bound} {percent} {count}".format(disc_bound=partition[1], percent=percent, count=this_stat)
			files[p][class_gp].write(outline+"\n")

for p in PRIMES:
	for class_gp in CLASS_GPS[p]:
		files[p][class_gp].close()

