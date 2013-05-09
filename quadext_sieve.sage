import itertools
import os, sys

if len(sys.argv) < 2:
	raise Exception("quadext_sieve.sage <upper bound>")

L.<zeta9> = NumberField(x^6 + x^3 + 1)

# This DOES include the end point.
def worthy_norm_sieve( upper_bound, L ):
	if not L.is_galois():
		raise NotImplementedError("The sieve requires L to be Galois.")
	upper_bound += 1 # so that we include the given upper bound
	sieve = bytearray( upper_bound >> 3 )
	# For simplicity, I use base-0 arrays.
	# Also '0' represents true while '1' represents false.
	# This is kinda stupid, but it's to save the pass over the bytearray
	# initializing everything to 1.
	# To access the ith bit (corresponding to integer value i+1), use sieve[i>>3] & (1 << (i & 7))
	# ie the sieve-index is (i >> 3) and the offset is (i & 7).
	for p in primes( upper_bound ):
		g_p = len( L(p).factor() )
		f_p = L.degree() / g_p
		# The below doesn't use for loops with range/xrange, since I want to use
		# Sage's integer type to avoid overflow issues.
		base = p
		while base < min(p^f_p, upper_bound):
			bad_integer = base
			while bad_integer < upper_bound:
				# sieve[bad_integer] = false
				index = bad_integer - 1
				index_high = index >> 3
				index_low  = index & 7
				sieve[index_high] |= (1 << index_low) # recall, 1 is false.
				bad_integer += p^f_p
			base += p
	good_integer = 1
	while good_integer < upper_bound:
		index = good_integer - 1
		index_high = index >> 3
		index_low  = index & 7
		if sieve[index_high] & (1 << index_low) == 0:
			# This means that good_integer is actually a good integer.
			yield good_integer
		#else: good_integer is not actually a good integer.
		good_integer += 1
	return

for norm in worthy_norm_sieve(UPPER_BOUND):
	print(norm)

