# Code to compute class groups given a list of number fields

import itertools
import os, sys
import time

if len(sys.argv) != 2:
	raise Exception("I need a list of files to process.")

FILE_LIST = None
with open(sys.argv[1]) as f:
	flines = f.readlines()
	FILE_LIST = map( lambda s: s.strip(), flines )

L.<zeta9> = NumberField(x^6 + x^3 + 1)

def parse_line( line ):
	els = line.split(':')
	norm = Integer(els[0])
	m_str = els[1]
	poly_tupled = map(tuple, eval(els[2]))
	poly = sum( map(lambda (c,p): c*(x^p), poly_tupled) )
	return norm, m_str, poly

@parallel
def pump_out_class_groups( infile_name, outfile_name ):
	infile  = open( infile_name,  'r' )
	outfile = open( outfile_name, 'w' )
	D = L.discriminant()
	D2 = D^2
	L4norm = L(4).norm()
	for line in infile:
		line = line.strip()
		# Ensure line isn't an error line
		if "ERROR" in line:
			outfile.write( "{oldstuff}:(No discriminant):(No class group)\n".format(oldstuff=line) )
			continue
		try:
			norm,m_str,polynomial = parse_line( line )
			K = NumberField( polynomial, 'sqrtm', cache=False )
			K_disc = K.discriminant()
			#if norm*D2 != K_disc and L4norm*norm*D2 != K_disc:
			if norm*D2 != K_disc:
				outfile.write( "{oldstuff}:ERROR norm did not match discriminant:(No class group)\n".format(oldstuff=line) )
				continue
			CG = K.class_group(proof=False) # NOTE THIS.
			elem_divisors = CG.elementary_divisors()
			elem_divisors_factored = map( lambda f: list(factor(f)), elem_divisors )
			elem_factors = []
			for l in elem_divisors_factored:
				elem_factors += l
			elem_divisors_sorted = sorted(elem_factors)
			prime_power_list = map( lambda (p,power): p^power, elem_divisors_sorted )
			outfile.write( "{oldstuff}:{disc}:{class_group}\n".format( oldstuff=line, disc=K_disc, class_group=prime_power_list ) )
		except Exception as e:
			outfile.write( "{oldstuff}:ERROR {msg}\n".format(oldstuff=line, msg=e) )
	outfile.close()
	infile.close()
	return

infiles = FILE_LIST
outfiles = map( lambda s: s.replace("polys","clsgps"), FILE_LIST )

total = len(FILE_LIST)
count = 0

start_time = time.time()

args_list = zip(infiles, outfiles)

for _ignore in pump_out_class_groups( args_list ):
	count += 1
	if count % 10 == 0:
		time_delta = time.time() - start_time
		perc = float(count) / total
		print("Finished {count}/{total} ({perc:.2f}%) T+{time_delta:.2f}s".format(count=count,total=total,perc=perc,time_delta=time_delta))

