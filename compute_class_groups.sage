# Code to compute class groups given a list of number fields

import itertools
import os

OUTDIR    = os.path.join( os.getcwd(), "Q_zeta9" )
INDIR     = os.path.join( os.getcwd(), "Q_zeta9" )
L.<zeta9> = NumberField(x^6 + x^3 + 1)
PART_SIZE = 10^4
START     = 10^0
END       = 10^8

# Test
PART_SIZE = 10^4
START     = 10^0
END       = 2*10^4

def parse_numberfield( string ):
	poly_coeffs_str = string.split(":")[2]
	poly_coeffs = map( Integer, str(filter(lambda c: c in "-0123456789,",poly_coeffs_str)).split(",") )
	poly = 0
	var = 1
	for c in poly_coeffs:
		poly += var*c
		var = var*x
	return NumberField( poly, 'z' )

@parallel
def pump_out_fields_in_range( dlow, dhigh ):
	infile_name  = os.path.join( INDIR,  "{dlow}-{dhigh}.polys.lst".format(dlow=dlow,dhigh=dhigh) )
	outfile_name = os.path.join( OUTDIR, "{dlow}-{dhigh}.clsgps.lst".format(dlow=dlow,dhigh=dhigh) )
	infile  = open( infile_name,  'r' )
	outfile = open( outfile_name, 'w' )
	for line in infile:
		line = line.strip()
		K = parse_numberfield( line )
		CG = K.class_group(proof=False) # NOTE THIS.
		elem_divisors = CG.elementary_divisors()
		elem_divisors_factored = flatten(map( factor, elem_divisors )) 
		elem_divisors_sorted = sorted(elem_divisors_factored)
		prime_power_list = map( lambda (p,power): p^power, elem_divisors_sorted )
		outfile.write( "{oldstuff}:{class_group}\n".format( oldstuff=line, class_group=prime_power_list ) )
	outfile.close()
	infile.close()
	return

def get_disc_partitions( low, high, size ):
	n_parts = (high - low + size - 1)/size # ceil((high - low)/size)
	for i in xrange( n_parts ):
		z = min( high, (i+1)*size )
		yield (low + i*size, z)

partitions = get_disc_partitions( START, END, PART_SIZE )
args_gen = itertools.imap( lambda (l,h): (l,h), partitions )

for _ignore1 in pump_out_fields_in_range( list(args_gen) ):
	print(_ignore1)

