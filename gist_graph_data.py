import itertools
import os
import sys

if len(sys.argv) != 7:
	raise Exception("gist_graph_data.py <primes_file> <data_file> <outdir> <disc_start> <disc_end> <disc_int>")

CLASS_GPS = {}
with open(sys.argv[1]) as f:
	CLASS_GPS = eval(f.read())
PRIMES = sorted(CLASS_GPS.keys())

DATA_FILE=open(sys.argv[2],'r')

OUTDIR    = sys.argv[3]
OUTFILE_T = os.path.join( OUTDIR, "graph_p={prime}_{class_gp}.dat" )

START     = Integer(sys.argv[4])
END       = Integer(sys.argv[5])
INTERVAL  = Integer(sys.argv[6])

def parse_line( line ):
	line = line.strip()
	# Ensure line isn't an error line
	if "ERROR" in line:
		return None
	try:
		# line is like "{norm}:{m}:{poly coefficients}:{discriminant}:{class group}"
		components = line.split(':')
		if len(components) < 5:
			return None
		discriminant = int(components[3])
		class_group_factors = eval(components[4])
		class_group_factors = tuple(class_group_factors)
		return discriminant, class_group_factors
	except:
		return None
	return None

def reduce_to_p_part( results, p ):
	new_results = {}
	for decomposition in results:
		ppart = filter( lambda n: n % p == 0, decomposition )
		ppart = tuple(ppart)
		new_results[ppart] = new_results.get(ppart,0) + results[decomposition]
	return new_results

def dump_results( results ):
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
threshold = START
for line in DATA_FILE:
	parsing = parse_line( line )
	if parsing is None:
		sys.stderr.write("Got invalid class group: {line}".format(line=line.strip()))
		continue
	this_disc, this_cls_gp = parsing
	if this_disc < START:
		continue
	if this_disc > threshold:
		dump_results( results )
		threshold += INTERVAL
	if this_disc > END:
		break
	results[this_cls_gp] = results.get(this_cls_gp,0) + 1

for p in PRIMES:
	for class_gp in CLASS_GPS[p]:
		files[p][class_gp].close()

