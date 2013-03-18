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
  , "prime_ramification" : A map of primes of ZZ which ramify in O_L to their factorization in O_L
  , "unit_adjusts"       : A collection of units of O_L which correspond
                           bijectively with the cosets of the 2-torsion of the
                           units of O_L.
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
  ramifying_primes = prime_factors(D)
  prime_ramification = {}
  for p in ramifying_primes:
    p_ideal_inQ = Q.fractional_ideal(p)
    p_ideal_inL = L.fractional_ideal(p)
    prime_ramification[p_ideal_inQ] = map( lambda (prime, power): prime, p_ideal_inL.factor() )
  U_L = L.unit_group()
  ngens = U_L.ngens()
  unit_adjust_vectors = [[]]
  for i in range(ngens):
    next_adjust_vectors  = map( lambda v: [0] + v, unit_adjust_vectors )
    next_adjust_vectors += map( lambda v: [1] + v, unit_adjust_vectors )
    unit_adjust_vectors = next_adjust_vectors
  unit_adjusts = map( U_L.exp, unit_adjust_vectors )
  return \
  { "discriminant"       : D
  , "prime_ramification" : prime_ramification
  , "unit_adjusts"       : unit_adjusts
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

def invert_norm( in_factors, L, prime_ramification ):
  """
  Let L/N be an extension of fields, wherein L and N are both number fields.

  For some ideal A of O_N, compute all ideals of O_L which have relative norm
  equal to A.
  (This might be extendable to any fractional ideals, but I don't care at this
  point.)

  The inputs are in_factors, which is the factorization of A, L, and
  prime_splittings, the information on how prime ideals split in L.

  The output is a list of prime factorizations of ideals of L which have relative
  norm A.

  NOTE: This has some optimizations for this program; in particular we're
  interested only in ideals which are square-free or <4>*square-free. This
  leads to some computational optimizations which are listed below and noted as
  they happen.
  """
  nonramifying_factors = []
  ramifying_factors = []
  itertools.ifilter( lambda t: not(t[0] in prime_ramification), in_factors )
  for prime, power in in_factors:
    if prime in prime_ramification:
      ramifying_factors.append( (prime,power) )
    else:
      nonramifying_factors.append( (prime,power) )
  #print( "Nonramifying factors: " + str(nonramifying_factors) )
  #print( "Ramifying factors   : " + str(ramifying_factors) )
  # The following for loop is an optimization for quadratic extensions.
  P2 = Q.fractional_ideal(2)
  for prime, power in nonramifying_factors:
    if power > 3 or (power > 1 and prime != Q.fractional_ideal(2)):
      # Everything that would be output by invert_norm will be divisible by
      # some ideal squared. So stop before doing all the hard work.
      return
  # The following for loop is an optimization for quadratic extensions.
  for prime, power in ramifying_factors:
    if power > len( prime_ramification[prime] ):
      # By the pigeon-hole principle, everything that would be output by
      # invert_norm must contain some upstairs prime which has power > 2,
      # implying the ideal is divisible by a square prime ideal.
      return
  # First, turn the nonramifying bits into ideals of M rather than N
  nonramified_factors = map( lambda (ideal, power): (L.fractional_ideal(ideal.gens()), power), nonramifying_factors )
  #print( "Nonramified factors : " + str(nonramified_factors) )
  # Second, ramify the ramifying bits.
  # The below uses partition_range_maxk with k=1 as an optimization for the
  # quadratic case; note that partition_range should be used in general.
  ramification_power_gens = itertools.imap( lambda (prime, power): partition_range_maxk(power, len(prime_ramification[prime]), 1), ramifying_factors )
  for ramification_powers in itertools.product( *ramification_power_gens ):
    # ramification_powers can be seen as a matrix wherein ramification_powers[i][j] is
    # the exponent on prime_ramification[ramifying_factors[i]][j]
    # That is, ramification_powers[i] is a list of powers for each of the
    # primes lying above p, where p = ramifying_factors[i][0].
    ramified_factors = []
    for i, ramification_powers_on_prime in enumerate( ramification_powers ):
      prime = ramifying_factors[i][0]
      ramified_factors += itertools.izip(prime_ramification[prime], ramification_powers_on_prime)
    #print( "Ramified factors    : " + str(ramified_factors) )
    yield Factorization( nonramified_factors + ramified_factors )

def expand_unit_adjusts( elt, unit_adjusts ):
  """
  Basically just map( lambda u: u*elt, unit_adjusts ).
  """
  for x in itertools.imap( lambda u: u*elt, unit_adjusts ):
    yield x

def generate_quadexts_withD( L, L_precomp, D ):
  """
  enumerate every quadratic extension of L (via L_precomp) with absolute
  discriminant equal to D exactly once.
  """
  if gcd( D, (L_precomp["discriminant"])^2 ) != (L_precomp["discriminant"])^2:
    return
  N = D/(L_precomp["discriminant"]^2)
  NI = Q.fractional_ideal(N)
  R.<x> = PolynomialRing(L)
  for relative_discriminant_factorization in invert_norm( list(NI.factor()), L, L_precomp["prime_ramification"] ):
    relative_discriminant_factorization = list(relative_discriminant_factorization)
    relative_discriminant = None
    if len(relative_discriminant_factorization) == 0:
      print("Got relative discriminant which is all of L")
      relative_discriminant = L.fractional_ideal(1)
    else:
      relative_discriminant = Factorization(relative_discriminant_factorization).expand()
    base_gen = relative_discriminant.gens_reduced()[0]
    for numbfld_gen in expand_unit_adjusts( base_gen, L_precomp["unit_adjusts"] ):
      if numbfld_gen == 1:
        continue
      yield L.extension( x^2 - numbfld_gen, 'm' )

L.<zeta9> = NumberField(x^6 + x^3 + 1)
precomps = precomputations(L)
D = 5*(precomps["discriminant"]^2)
for k in precomps:
  print(str(k) + " => \t" + str(precomps[k]))
print

Ks = []
for K in generate_quadexts_withD( L, precomps, 5*(precomps["discriminant"]^2) ):
  Ks.append(K)

print("Got {n} number fields:".format(n=len(Ks)))
for K in Ks:
  print("\t{K}".format(K=K))
print

print("Discriminants. Expect D={D}".format(D=D))
for K in Ks:
  print("{K} -> disc(K)={D}".format(K=K,D=K.absolute_field('a').discriminant()))

bad = 0
for K1,K2 in itertools.combinations(Ks,2):
  if K1.is_isomorphic(K2):
    print("ISOMORPHISM:")
    print("\t" + str(K1))
    print("\t" + str(K2))
    bad += 1
print("# pairs which are isomorphic: {bad}".format(bad=bad))

