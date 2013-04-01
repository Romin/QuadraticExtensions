# Code to enumerate all quadratic extensions of a given number field
# (generally referred to as L)
#
# One requirement is that the field's class number must be 1 for this to work.

import itertools

# For my own sanity.
Q = NumberField(x,'a')

def precomputations( L, verify=True ):
	"""
	We need to know a few things about L that we can precompute before
	enumerating fields. This function does those computations, and returns a
	dictionary keyed by strings with those values, as follows:
	{ "discriminant"       : Absolute discriminant of L
	, "is_galois"          : Boolean that is true iff L is galois
	, "unit_adjusts"       : A collection of units of O_L which correspond
													 bijectively with the cosets of the 2-torsion of the
													 units of O_L.
	, "mod4"               : The ring O_L / 4O_L
	, "squares_mod4"       : The set of elements x^2 where x is in O_L / 4O_L.
	}
	
	Also, if verify is True, then L is verified to have class number 1. If it
	does not, None is returned instead of the above dictionary.
	If verify is false, then L is *not* verified to have class number 1, and the
	results of this are undefined when L does not have class number 1.
	"""
	if verify:
		if L.class_number() != 1:
			return None
	D = L.discriminant()
	is_galois = L.is_galois()
	U_L = L.unit_group()
	ngens = U_L.ngens()
	unit_adjust_vectors = [[]]
	for i in range(ngens):
		next_adjust_vectors  = map( lambda v: [0] + v, unit_adjust_vectors )
		next_adjust_vectors += map( lambda v: [1] + v, unit_adjust_vectors )
		unit_adjust_vectors = next_adjust_vectors
	unit_adjusts = map( U_L.exp, unit_adjust_vectors )
	# TODO: This method of creating a set of all squares mod 4 is pretty stupid
	# and could probably be improved.
	O_L = L.ring_of_integers()
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
	mod4_elts = itertools.imap( lambda elts: sum(map(prod,elts)), superpaired )
	mod4_squares = Set( [O_L_4O_L(a*a) for a in mod4_elts] )
	return \
	{ "discriminant"       : D
	, "is_galois"          : is_galois
	, "unit_adjusts"       : unit_adjusts
	, "mod4"               : O_L_4O_L
	, "squares_mod4"       : mod4_squares
	}

def partition_range( S, n ):
	"""
	Uniquely generate all lists `lst' with the property that sum(lst)=S and
	len(lst) = n and x in S implies x >= 0.

	I assume n > 0 and S >= 0.
	"""
	if n <= 0 or S > 0:
		raise Exception("partition_range got bad arguments.")
	if n == 1:
		yield [S]
	else:
		for i in xrange(S+1):
			prefix = [i]
			for suffix in partition_range(S-i, n-1):
				yield prefix + suffix

def partition_range_maxk( S, n, k ):
	"""
	A variant of partition_range which has useful optimizations for this program.

	It uniquely generates all lists `lst' with the property that sum(lst)=S and
	len(lst) = n and x in S implies x >= 0, as before, with the additional
	constraint that x in S implies x <= k.

	Note that this *ASSUMES* that the given S, n, k combination is feasible.
	In particular, 0 <= S <= n*k, n > 0, and k > 0.
	"""
	if S < 0 or S > n*k or n <= 0 or k <= 0:
		raise Exception("partition_range_maxk got bad arguments.")
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
			for suffix in partition_range_maxk( S-i, n-1, k ):
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
	"""
	if len( in_factors ) == 0:
		yield Factorization([])
		return
	downstairs_ideal, downstairs_power = in_factors[0]
	upstairs_factorization = L.fractional_ideal(downstairs_ideal.gens()).factor()
	upstairs_ideals = map( lambda (ideal, power): ideal, upstairs_factorization )
	if L_precomp["is_galois"]:
		residue_class_degree = upstairs_ideals[0].residue_class_degree()
		if not residue_class_degree.divides( downstairs_power ):
			return
		# The 1 here is the `k' argument to partition_range_maxk
		if downstairs_power/residue_class_degree > len(upstairs_ideals)*1:
			# We'd pass invalid arguments to partition_range_maxk.
			# Plus the pigeon-hole principle guarantees that we cannot have a
			# square-free ideal.
			return
		for other_factors in invert_norm( in_factors[1:], L, L_precomp ):
			for exponents in partition_range_maxk( downstairs_power/residue_class_degree, len(upstairs_ideals), 1 ):
				this_factor = Factorization( zip(upstairs_ideals, exponents) )
				yield this_factor * other_factors
	else:
		raise NotImplementedError("DON'T GIVE NONGALOIS FIELDS YET")
	return

def expand_unit_adjusts( elt, unit_adjusts ):
	"""
	Basically just map( lambda u: u*elt, unit_adjusts ).
	"""
	for x in itertools.imap( lambda u: u*elt, unit_adjusts ):
		yield x

def generate_quadexts_withD( L, L_precomp, D ):
	"""
	enumerate every quadratic extension of L (via L_precomp) with absolute value
	of the absolute discriminant equal to D exactly once.
	"""
	if not ((L_precomp["discriminant"])^2).divides( D ):
		return
	N = D/(L_precomp["discriminant"]^2)
	NI = Q.fractional_ideal(N)
	R.<x> = PolynomialRing(L)
	relative_discriminant_factorization_generator = invert_norm( list(NI.factor()), L, L_precomp )
	# The flag here being True indicates that the ideals in the associated
	# generator will have the 4 pulled out of them.
	flag_and_gens = [ (False, relative_discriminant_factorization_generator) ]
	if L(4).norm().divides( N ):
		NI4 = Q.fractional_ideal( N / (L(4).norm()) )
		relative_discriminant_factorization_generator = invert_norm( list(NI4.factor()), L, L_precomp )
		flag_and_gens.append( (True, relative_discriminant_factorization_generator) )
	for flag, relative_discriminant_factorization_generator in flag_and_gens:
		for relative_discriminant_factorization in relative_discriminant_factorization_generator:
			relative_discriminant_factorization = list(relative_discriminant_factorization)
			relative_discriminant = None
			if len(relative_discriminant_factorization) == 0:
				print("Got relative discriminant which is all of L")
				relative_discriminant = L.fractional_ideal(1)
			else:
				relative_discriminant = Factorization(relative_discriminant_factorization).expand()
			base_gen = relative_discriminant.gens_reduced()[0]
			for numbfld_gen in expand_unit_adjusts( base_gen, L_precomp["unit_adjusts"] ):
				if numbfld_gen == 1: # Skip 1 because it's dumb
					continue
				numbfld_gen_mod4 = L_precomp["mod4"]( numbfld_gen )
				if flag: # We want numbfld_gen _IS_NOT_ a square mod 4
					if numbfld_gen_mod4 in L_precomp["squares_mod4"]:
						print("D={D}, N={N}, {x} is a square mod 4".format(x=numbfld_gen, D=D, N=N/L(4).norm()))
						continue
				else: # We want numbfld_gen _IS_ a square mod 4
					if numbfld_gen_mod4 not in L_precomp["squares_mod4"]:
						print("D={D}, N={N}, {x} is not a square mod 4".format(x=numbfld_gen, D=D, N=N))
						continue
				yield L.extension( x^2 - numbfld_gen, 'm' )
	return

L.<zeta9> = NumberField(x^6 + x^3 + 1)
precomps = precomputations(L)
#D = 3*(L(4).norm())*(precomps["discriminant"]^2)
D = 3*(precomps["discriminant"]^2)
print("Pre-computations completed. Continuing with enumeration.")
print

Ks = []
for K in generate_quadexts_withD( L, precomps, D ):
	Ks.append(K)

print("Got {n} number fields:".format(n=len(Ks)))
for K in Ks:
	print("\t{K}".format(K=K))
print

print("Discriminants. Expect D={D}".format(D=D))
for K in Ks:
	print("{K} -> disc(K)={D}".format(K=K,D=K.absolute_field('a').discriminant()))
print

bad = 0
for K1,K2 in itertools.combinations(Ks,2):
	if K1.is_isomorphic(K2):
		print("ISOMORPHISM:")
		print("\t" + str(K1))
		print("\t" + str(K2))
		bad += 1
print("# pairs which are isomorphic: {bad}".format(bad=bad))

