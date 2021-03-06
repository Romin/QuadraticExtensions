# Code to enumerate all quadratic extensions of a given number field
# (generally referred to as L)
#
# One requirement is that the field's class number must be 1 for this to work.
#
# Also, though not a theoretical requirement for this algorithm, there are
# optimizations that exist when L is Galois, and only those optimizations are
# implemented; the non-Galois case is not yet implemented.

import itertools
import os, sys
import time

# Sage has QQ, but that doesn't have NumberField capabilities, so this serves
# as our canonical field of rationals.
Q = NumberField(x,'a')
R.<y> = PolynomialRing(Q)

# The parameters of this computation
OUTDIR    = os.path.join( os.getcwd(), "Q_zeta9" )
L.<zeta9> = NumberField(x^6 + x^3 + 1)
PART_SIZE = 10^4
START     = 1
END       = 10^9

NCPUS     = None # None is the default

# Globals that depend on the parameters
S.<z>     = PolynomialRing(L)

def precomputations( L, **kwargs ):
	"""
	We need to know a few things about L that we can precompute before
	enumerating fields. This function does those computations, and returns a
	dictionary keyed by strings with those values, as follows:
	{ "discriminant"          : Absolute discriminant of L
	, "is_galois"             : Boolean that is true iff L is galois
	, "unit_adjusts"          : A collection of units of O_L which correspond
	                            bijectively with the cosets of the 2-torsion of the
	                            units of O_L.
	, "mod4_elts"             : Representatives of every coset of O_L / 4O_L
	, "unit_to_squaremod4"    : A function mapping x in O_L to u, a unit of O_L,
	                            so that x*u is a square mod 4, or None if no such unit
	                            exists.
	, "4norm"                 : L(4).norm()
	, "prime_splits"          : A map of primes of Q to information on their factoring in L.
	                            p -> [(q_i, e_i, f_i)] where p = prod(q_i^e_i) and
	                            q_i.norm() = p^f_i
															upper_bound must be specified for this
															precomputation to occur. p <= upper_bound are
															computed.
	, "norms"                 : A list of integers up to upper_bound from
															norms_file which are norms of ideals of O_L.
															Both upper_bound and norms_file must be specified
															for this precomputation to occur.
	}

	Keyword parameters:

	upper_bound is a positive integer which is an upper bound on lists of numbers
	which are computed.
	If upper_bound is None, then none of the precomputations which require
	upper_bound will occur.

	norms_file is a string which is a filename which refers to a file containing
	a list, one element per line, of positive integers which are absolute values
	of norms of ideals of L.
	If norms_file is None, then precomputations relying on norms_file will not
	occur.
	
	verify is a boolean which, when True, includes a precomputation in which L is
	verified to have class number 1. If L does not have class number 1, an
	exception is raised.
	If verify is false, then L is *not* verified to have class number 1, and the
	results of this are undefined when L does not have class number 1.
	"""
	if kwargs.get("verify",True):
		print("\tVerifying L has class number 1 ...")
		if L.class_number() != 1:
			raise Exception("L does not have class number 1.")
	print("\tComputing L's discriminant...")
	D = L.discriminant()
	print("\tComputing if L is Galois...")
	is_galois = L.is_galois()
	print("\tComputing unit adjusts...")
	U_L = L.unit_group()
	O_L = L.ring_of_integers()
	ngens = U_L.ngens()
	unit_adjust_vectors = [[]]
	for i in range(ngens):
		next_adjust_vectors  = map( lambda v: [0] + v, unit_adjust_vectors )
		next_adjust_vectors += map( lambda v: [1] + v, unit_adjust_vectors )
		unit_adjust_vectors = next_adjust_vectors
	unit_adjusts = map( U_L.exp, unit_adjust_vectors )
	unit_adjusts = [ O_L(u) for u in unit_adjusts ]
	print("\tComputing squares mod 4O_L...")
	# TODO: This method of creating a set of all squares mod 4 could probably be
	# improved.
	# For example, with L=Q(zeta9), mod4_elt_vectors gets to size 4096, but
	# mod4_squares has size only 64.
	# It seems to be the case that a lifting of elements of O_L/2O_L to O_L/4O_L
	# will give us precisely the squares mod 4O_L.
	O_L_4O_L = O_L.quotient( O_L.ideal(4), 'q' )
	basis = O_L.basis()
	mod4_elt_vectors = [[]]
	for i in range(len(basis)):
		next_elt_vectors  = []
		for j in range(4):
			next_elt_vectors += map( lambda v: [j] + v, mod4_elt_vectors )
		mod4_elt_vectors = next_elt_vectors
	paired = itertools.izip( mod4_elt_vectors, itertools.repeat(basis) )
	superpaired = itertools.starmap( zip, paired )
	mod4_elts = list(itertools.imap(lambda elts: sum(map(prod,elts)),superpaired))
	mod4_squares = Set( [O_L_4O_L(a*a) for a in mod4_elts] )
	print("\tComputing unit-adjust-to-square-mod-4 hash table...")
	def idx( x ):
		coords = x.list()
		coords = reversed(coords)
		coords = [ c % 4 for c in coords ]
		index = 0
		for c in coords:
			index *= 4
			index += c
		return index
	unit_to_squaremod4_hashtable = [None for elt in mod4_elts]
	for elt in mod4_elts:
		for u in unit_adjusts:
			if O_L_4O_L(elt*u) in mod4_squares:
				unit_to_squaremod4_hashtable[idx(elt)] = u
	def unit_to_squaremod4( x ): # Closures are wonderful, wonderful things.
		return unit_to_squaremod4_hashtable[idx(x)]
	print("\tComputing 4's norm...")
	norm4 = L(4).norm()
	prime_splits = {}
	if kwargs.get("upper_bound",None) is not None:
		print("\tComputing prime splitting information...")
		upper_bound = ceil( sqrt(kwargs["upper_bound"]) ) # Larger primes are unlikely to appear often
		primes_list = prime_range(upper_bound+1)
		primes_list_list = []
		part_size=10^3
		n_parts = ceil(len(primes_list)/part_size)
		for i in range( n_parts ):
			if i == n_parts-1:
				primes_list_list.append( primes_list[i*part_size:] )
			else:
				primes_list_list.append( primes_list[i*part_size:(i+1)*part_size] )
		primes_list = None
		@parallel(ncpus=NCPUS)
		def gen_prime_info( primes_list ):
			prime_splits = {}
			for p in primes_list:
				up = L(p)
				up_factored = up.factor()
				if is_galois:
					residue_class_degree = L.ideal(up_factored[0][0]).residue_class_degree()
					f = lambda q: residue_class_degree # constant function
				else:
					f = lambda q: L.ideal(q).residue_class_degree()
				prime_splits[p] = [(q_i, e_i, f(q_i)) for (q_i, e_i) in up_factored]
			return prime_splits
		for args, result in gen_prime_info(primes_list_list):
			for k in result:
				prime_splits[k] = result[k]
			result = None
	norms = []
	if kwargs.get("upper_bound",None) is not None and kwargs.get("norms_file",None) is not None:
		upper_bound = kwargs["upper_bound"]
		norms_filename = kwargs["norms_file"]
		print("\tReading norms from \"{norms_file}\"".format(norms_file=norms_filename))
		norms_file = open( norms_filename, 'r' )
		for line in norms_file:
			try:
				norm = Integer(line.strip()) # Pray there are no leading zeroes.
			except:
				continue
			if norm > upper_bound:
				break
			norms.append( norm )
		norms_file.close()
	return \
	{ "discriminant"       : D
	, "is_galois"          : is_galois
	, "unit_adjusts"       : unit_adjusts
	, "mod4_elts"          : mod4_elts
	, "unit_to_squaremod4" : unit_to_squaremod4
	, "4norm"              : norm4
	, "prime_splits"       : prime_splits
	, "norms"              : norms
	}

def integer_partition_range_maxk( S, n, k ):
	"""
	A variant of partition_range which has useful optimizations for this program.

	It uniquely generates all lists `lst' with the property that sum(lst)=S and
	len(lst) = n and x in S implies x >= 0, as before, with the additional
	constraint that x in S implies x <= k.

	Note that this *ASSUMES* that the given S, n, k combination is feasible.
	In particular, 0 <= S <= n*k, n > 0, and k > 0.
	"""
	if S < 0 or S > n*k or n <= 0 or k < 0:
		return # There are no such `lst'
	if n == 1:
		yield [S]
	else:
		# If I put i into this container with i < S - (n-1)*k, then
		# I'll have S - i > S - S + (n-1)*k = (n-1)*k things to put into n-1 boxes,
		# which contradicts the assumption I make.
		# If I put S - (n-1)*k into this container, then
		# I'll have S - S + (n-1)*k = (n-1)*k things to put into n-1 boxes, which
		# is valid.
		minm_i = max( S - (n-1)*k, 0 )
		for i in xrange( minm_i, min(S,k)+1 ):
			prefix = [i]
			for suffix in integer_partition_range_maxk( S-i, n-1, k ):
				yield prefix + suffix

def invert_norm( in_factors, L, L_precomp ):
	"""
	Let L/N be an extension of fields, wherein L and N are both number fields.

	For some ideal A of O_N, compute all ideals of O_L which have relative norm
	equal to A.
	(This might be extendable to any fractional ideals in the sense that some
	factors have negative exponents, but I don't care at this point.)

	The inputs are in_factors, which is the factorization of A, L, and
	prime_splittings, the information on how prime ideals split in L.

	The output is a list of prime factorizations of ideals of L which have relative
	norm A.

	NOTE: This has some optimizations for this program; in particular we're
	interested only in ideals which are square-free. This leads to some
	computational optimizations which are listed below and noted as they happen.
	- partition_range_maxk is used with k=1 to ensure squarefree-ness instead of
		partition_range.
	- the non-galois case is not implemented, since we're only interested in
		Q(z_9), which is galois.
	- ideal generators are used in lieu of actual ideals, since we're only
		interesting in L which have trivial class group/are PIDs.
	"""
	if len( in_factors ) == 0:
		yield Factorization([])
		return
	downstairs_ideal_gen, downstairs_power = in_factors[0]
	upstairs_factorization_info = L_precomp.get("prime_splits",{}).get(downstairs_ideal_gen,None)
	if upstairs_factorization_info is None:
		upstairs_factorization = L(downstairs_ideal_gen).factor()
		if L_precomp["is_galois"]:
			f_p = L.ideal(upstairs_factorization[0][0]).residue_class_degree()
			upstairs_factorization_info = [ (q,e,f_p) for (q,e) in upstairs_factorization ]
		else:
			f = lambda q: L.ideal(q).residue_class_degree()
			upstairs_factorization_info = [ (q,e,f(q)) for (q,e) in upstairs_factorization ]
	upstairs_ideal_gens = [ q for (q,e,f) in upstairs_factorization_info ]
	if L_precomp["is_galois"]:
		residue_class_degree = upstairs_factorization_info[0][2]
		if not residue_class_degree.divides( downstairs_power ):
			return
		# The 1 here is the `k' argument to partition_range_maxk
		if downstairs_power/residue_class_degree > len(upstairs_ideal_gens)*1:
			# We'd pass invalid arguments to partition_range_maxk.
			# Plus the pigeon-hole principle guarantees that we cannot have a
			# square-free ideal.
			return
		for other_factors in invert_norm( in_factors[1:], L, L_precomp ):
			for exponents in integer_partition_range_maxk( downstairs_power/residue_class_degree, len(upstairs_ideal_gens), 1 ):
				this_factor = Factorization( zip(upstairs_ideal_gens, exponents) )
				yield this_factor * other_factors
	else:
		raise NotImplementedError("Currently only L which are Galois are supported by invert_norm.")
	return

def generate_quadexts_with_norm( L, L_precomp, norm ):
	"""
	enumerate every quadratic extension of L (via L_precomp) with absolute value
	of the norm of the relative discriminant over L equal to norm exactly once.
	"""
	N = norm
	N_factored = N.factor()
	danger_zone = False
	if norm > L_precomp["upper_bound"] / L_precomp["4norm"]:
		danger_zone = True
	if danger_zone:
		if norm % 4 != 1:
			# There won't be anything that we're interested in.
			return
	for relative_discriminant_gen_factorization in invert_norm( list(N_factored), L, L_precomp ):
		base_gen = None
		if len(relative_discriminant_gen_factorization) == 0:
			# Pick 1 as the generator to simplify the filtering of square units
			# later. (1 will be the only potential square unit.)
			base_gen = L(1)
		else:
			base_gen = relative_discriminant_gen_factorization.expand()
		unit_adjusts = []
		if danger_zone:
			unit_adjust = L_precomp["unit_to_squaremod4"]( base_gen )
			if unit_adjust is not None:
				unit_adjusts = [unit_adjust]
		else:
			unit_adjusts = L_precomp["unit_adjusts"]
		for unit in unit_adjusts:
			numbfld_gen = unit*base_gen
			if numbfld_gen == L(1):
				continue
			abs_poly = numbfld_gen.minpoly()(x^2)
			yield numbfld_gen, abs_poly
	return

@parallel(ncpus=NCPUS)
def pump_out_fields( bounds, norms, L, L_precomps ):
	bounds_str = "{low}-{high}".format(low=bounds[0], high=bounds[1])
	outfile = open( "{partition_id}.polys.lst".format(partition_id=bounds_str), 'w' )
	for norm in norms:
		D = norm * precomps["discriminant"]^2
		try:
			for m, abs_poly in generate_quadexts_with_norm( L, precomps, norm ):
				line = "{norm}:{m}:{coefficients}".format( norm=norm, m=m, coefficients=abs_poly.coeffs() )
				outfile.write( line+"\n" )
		except Exception as e:
			outfile.write( "{norm}:ERROR:\"{msg}\"\n".format(norm=norm,msg=e) )
	outfile.close()
	return

def get_partitions( low, high, size ):
	n_parts = ceil( (high - low)/size )
	for i in xrange( n_parts ):
		z = min( high, low + (i+1)*size - 1 )
		yield (low + i*size, z)

def get_partitioned_norms( norms, partition_bounds ):
	partitions = {}
	for bound in partition_bounds:
		partitions[bound] = []
	i = 0
	j = 0
	while i < len(norms) and j < len(partition_bounds):
		norm = norms[i]
		bound = partition_bounds[j]
		if norm > bound[1]:
			j += 1
		elif norm < bound[0]:
			i += 1
		else:
			partitions[bound].append(norm)
			i += 1
	return [partitions[bound] for bound in partition_bounds]

os.chdir( OUTDIR )

print("Beginning precomputations on L...")
precomps = precomputations(L, upper_bound=END, norms_file="norms.lst")

precomps["upper_bound"] = END

partitions = list(get_partitions(START, END, PART_SIZE))
partitioned_norms = get_partitioned_norms( precomps["norms"], partitions )

arg_list = [(bounds, norms_part, L, precomps) for (bounds, norms_part) in zip(partitions, partitioned_norms)]

START_TIME = time.time()
print("Generating number fields... (START TIME={start_time})".format(start_time=START_TIME))
for args, _ignore in pump_out_fields( arg_list ):
	args = args[0] # don't want keyword args
	bounds = args[0]
	norms = args[1]
	time_delta = time.time() - START_TIME
	print("Finished range [{low},{high}] ({n_norms} norms) @T+{time}s".format(low=bounds[0], high=bounds[1], n_norms=len(norms), time=time_delta))

