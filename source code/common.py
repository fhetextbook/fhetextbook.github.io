import numpy as np
from numpy.polynomial import Polynomial
from fastcore.foundation import patch_to
from os import sys
import random
import math
from sympy import nextprime, mod_inverse, primitive_root
import copy
import secrets 

np.set_printoptions(suppress=True)



def poly_rescale_P_over_Q(raw_poly: Polynomial, P: int, Q: int, N: int) -> Polynomial:
	"""
	rescale: round((P/Q) * coeff) coefficient-wise, then reduce mod Q.
	Output is in R_Q.
	"""
	raw = [int(x) for x in np.array(raw_poly.coef, dtype=object).reshape(-1)]
	scaled = [round_div(P * c, Q) for c in raw]	 # big-int rounding
	return poly_mod(Polynomial(np.array(scaled, dtype=object)), N, Q)		
	

def rand_odd(lo: int, hi: int) -> int:
	n = lo + secrets.randbelow(hi - lo)
	if n % 2 == 0:
		n += 1
	if n >= hi:
		n -= 2
	return n


def is_prime(n, k=5):
	"""
	Miller-Rabin test with a static cache to store results.
	"""
	# 1. Initialize cache if it doesn't exist (Static Storage)
	if not hasattr(is_prime, "cache"):
		is_prime.cache = {}

	# 2. Check Cache
	if n in is_prime.cache:
		return is_prime.cache[n]

	# --- Standard Checks ---
	if n < 2: 
		is_prime.cache[n] = False
		return False
	if n == 2 or n == 3: 
		is_prime.cache[n] = True
		return True
	if n % 2 == 0: 
		is_prime.cache[n] = False
		return False

	# 3. Small Prime Filter (Optimization)
	for p in [3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]:
		if n == p: 
			is_prime.cache[n] = True
			return True
		if n % p == 0: 
			is_prime.cache[n] = False
			return False

	# --- Miller-Rabin Logic ---
	r, d = 0, n - 1
	while d % 2 == 0:
		r += 1
		d //= 2

	# Deterministic bases for 64-bit range, random otherwise
	if n < 18446744073709551616:
		bases = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]
	else:
		bases = [random.randrange(2, n - 1) for _ in range(k)]

	for a in bases:
		if a >= n: break
		x = pow(a, d, n)
		if x == 1 or x == n - 1:
			continue
		for _ in range(r - 1):
			x = pow(x, 2, n)
			if x == n - 1:
				break
		else:
			# Composite found
			is_prime.cache[n] = False
			return False
			
	# Probably Prime
	is_prime.cache[n] = True
	return True



def pick_primes(bits: int = 30, k: int = 5, max_span: int = 1 << 20) -> list[int]:
	lo, hi = 1 << (bits - 1), 1 << bits  # exactly 'bits' bits
	while True:
		center = rand_odd(lo, hi)
		w_lo = max(lo, center - max_span // 2)
		w_hi = min(hi, center + max_span // 2)
		if w_hi - w_lo < 10:  # degenerate window; retry
			continue

		primes = set()
		for _ in range(200_000):
			cand = rand_odd(w_lo, w_hi)
			if is_prime(cand):
				primes.add(cand)
				if len(primes) == k:
					ps = sorted(primes)
					if ps[-1] - ps[0] <= max_span:
						return ps


def pick_ntt_primes(bits: int = 30, num: int = 5, N: int = 4, p_mod: int = None, max_span: int = 1 << 20) -> list[int]:
    """
    Picks 'k' primes such that:
      1. q = 1 mod 2N (NTT requirement)
      2. (Optional) q = 1 mod P (if P is not None)
      3. Primes are within 'max_span' of the largest picked prime.
    """
    ntt_mod = 2 * N
    
    if p_mod is not None:
        # Enforce both 1 mod 2N and 1 mod P -> LCM(2N, P)
        M = (ntt_mod * p_mod) // math.gcd(ntt_mod, p_mod)
    else:
        # Only enforce 1 mod 2N
        M = ntt_mod
    
    primes = []
    
    # Start looking near 2^bits
    start = (1 << bits)
    
    # Ensure start = 1 mod M
    start -= (start % M)
    start += 1
    
    candidate = start
    
    while len(primes) < num:
        # Search downwards to stay close to 2^bits
        candidate -= M
        
        # 1. Bit-length safety check
        if candidate < (1 << (bits - 1)):
            raise ValueError(f"Could not find {num} NTT primes in this bit range satisfying requirements.")
            
        # 2. Max-Span check
        # Since we search downwards, primes[0] is the largest.
        if len(primes) > 0 and (primes[0] - candidate > max_span):
            raise ValueError(f"Could not find {num} primes within max_span {max_span}")
            
        # 3. Fast Primality Check (Miller-Rabin)
        if is_prime(candidate):
            primes.append(candidate)
            
    return sorted(primes)
	

def centered_mod(a, q):
	a %= q
	half = q // 2
	# even q: [-q/2, q/2), odd q: [-(q//2), q//2]
	if (q % 2 == 0 and a >= half) or (q % 2 == 1 and a > half):
		a -= q
	return a



def centered_mod_arr(a, q):
	a = np.array(a, dtype=object)
	a = a % q
	half = q // 2
	# even q: [-q/2, q/2) / odd q: [-(q//2), q//2]
	if q % 2 == 0:
		return np.where(a >= half, a - q, a)
	else:
		return np.where(a > half, a - q, a)



def round_div(x, d):
	"""Nearest integer to x/d, ties away from 0 (integer-only)."""
	#x = int(x); d = int(d)
	if d <= 0:
		raise ValueError("scale must be positive")
	if x >= 0:
		return (x + d // 2) // d
	return -((-x + d // 2) // d)



def is_power_of_2(N):
	n_checker = N
	while n_checker > 1:
		if n_checker % 2 == 1:
			return False
		n_checker /= 2 
	return True



def egcd(a, b):
	if a == 0:
		return (b, 0, 1)
	else:
		g, y, x = egcd(b % a, a)
		return (g, x - (b // a) * y, y)



def mod_inv(a, m):
	g, x, y = egcd(a, m)
	if g != 1:
		#raise Exception('modular inverse does not exist')
		return -1
	else:
		return x % m




def init_P(N):

	M = N*2
	if not is_power_of_2(N):
		print('[Error] N is not a power of 2: ' + str(N))
		sys.exit()

	# Fermit's Little theorem: if A and p are co-primes, then A^{p-1} ≡ 1 mod p 
	# ==> If p is a prime, then any number A in [0, p - 1] has co-prime with p, thus any number A in [0, p - 1] has a multiplicative inverse, which is A^{p-2} mod p, and A * A^{p-2} = A^{p-1} ≡ 1 mod p

	# Find a modulus P that guarnatees the existence of a primitive M-th root of unity and the inverse of N
	# If P is a prime, then the inverse of N always exists
	# A primitive (M=2N)-th root of unity exists if M=2N divides P - 1 ???
	#  X is the primitive M-th root of unity if X^M ≡ 1 mod P AND Order(X) = M (where M = 2N)

	P = nextprime(M)  # Find a prime P bigger than M(=2N)
	# We find P such that there exists the primitive M-th root of unity in modulo P. According to Fermat's Little Theorem (A^{P-1} ≡ 1 mod P), each element in P has an order that divides P-1 (i.e., the multiplicative modulo group size), and all elements in the group collectively cover all divisors of P-1. Among these elements, given Ord(A) = M for some element A, that's the primitive M-th root of unity modulo P (as the definition of primitive M-th root of unity: any M-th root of unity C, solution for X^M ≡ 1, such that Ord(C) = M). Therefore, if P-1 (the multiplicative modulo group size) is some multiple of M, then A^{P-1} = A^{M*k} = (A^k)^M ≡ 1 mod P, then some element B=A^k in this group will have Ord(B) = M. As a conclusion, if (P - 1) mod M == 0 (i.e., M divides the size of the multiplicative group of P), then some elements in this group will have an order M, which will be the primitive M-th root of unity. 

	while (P - 1) % M != 0:
		P = nextprime(P)
	return P


def print_cyclotomic_polynomial_info(N=None, P=None, is_complex=False):
	M = N*2
	if  is_complex:
		N_inverse = 1/N
		root_values = []
		primitive_root_of_unity = np.exp(2 * np.pi * 1j / M) 
	else:
		N_inverse = mod_inv(N, P) # If P is a prime, the multiplicative inverse of N (< P) is guaranteed to exist (by Fermit's Little theorem)
		mult_group_generator = primitive_root(P)  # prim_root is a number such that Order(prim_root) = totient(p) mod p
		primitive_root_of_unity = pow(mult_group_generator, (P - 1) // M, P) # The (M=2N)-th primitive root of unity (i.e., solution for both X^N + 1 = 0 mod p AND X^{2N} = 1 mod P)

	root_values = []

	# [OPTION 1]: A brute-force way of finding the primitiev (M=2N)-th root of unity (i.e., the solution for X^N + 1 = 0 mod P)
	'''
	# After then, we find X ∈ {1,2,...P-1} whose order is 2N. That value is the primitive P-th root of unity.
	while len(root_values) == 0:

		# Find P such that X^{P-1} ≡ 1 mod P  (Fermat's little theorem) AND M=2N divides P-1, which ensures that the order of each X ∈ {1,2,...P-1} mod P is either 2N or some factor of 2N

		N_inverse = ModInv(N, P) # If P is a prime, the multiplicative inverse of N (< P) is guaranteed to exist (by Fermit's Little theorem)
		if not is_random_param and N_inverse <= 0:
			print("[Error] N=" + str(N) + "has no inverse modulo P=" + str(P))
			sys.exit()
		primitive_root_of_unity = -1
		for i in range(P):
			result = False
			if is_random_param:
				result = EvaluateRandomCyclotomicPolynomial(i) % P == 0 # If x^N ≡ -1 mod p (i.e., X^N + 1 ≡ 0 mod p)
			else:
				result = EvaluateTargetedCyclotomicPolynomial(i) % P == 0
				print  
			if result: # If i is the solution for X^N + 1 (i.e., M/2-th cyclotomic polynomial) 
				if (i ** (M)) % P == 1 and (M % 2 == 1 or (i ** N) % P != 1): # If Ord(i) = M
					root_values.append(i)	# Then, i is the primitive M-th root of unity
		if len(root_values) == 0 and not is_random_param: 
			print("Modulus P=" + str(P) + " has no primitive " + str(M) + "-th root of unity")
			sys.exit(0)

	primitive_root_of_unity = min(root_values)
	'''

#  '''
	### [OPTION 2]: An alternative efficient way of computing the primitive (M=2N)-th root of unity
	next_root = primitive_root_of_unity
	if is_complex:
		for i in range(N):
			root_values.append(next_root)
			next_root = next_root * (primitive_root_of_unity ** 2)
	else:
		while True:
			root_values.append(next_root)
			next_root = next_root * (primitive_root_of_unity ** 2) % P
			if next_root == primitive_root_of_unity:
				break
#  '''

	print("==================== PARAMETERS ====================")
	print("N: " + str(N))
	if not is_complex:
		print("Plaintext modulus: " + str(P))
	#print("Scale Factor: " + str(scale))
	print("Selected primitive " + str(M) + "-th root of unity: " + str(primitive_root_of_unity))
	print("- " + str(N) + " roots of the " + str(M) + '-th cyclotomic polynomial : ' + str(root_values))
	print("- The inverse of N=" + str(N) + " is " + str(N_inverse))
	print()
	for root in root_values:
		generated_roots = []
		coprimes = []
		for exp in range(M):
			if exp != 0 and math.gcd(exp, M) == 1:
				coprimes.append(exp)
				if is_complex:
					generated_roots.append((int(root) ** exp))
				else:
					generated_roots.append((int(root) ** exp) % P)
		print('- Roots of the ' + str(M) + '-th cyclotomic polynomial generated by the co-prime-' + str(coprimes) + '-exponented primitive ' + str(M) + '-th root of unity ' + str(root) + ' : ' + str(generated_roots))
	print()
	for root in root_values:
		generated_roots = []
		exponents = []
		for exp in range(M):
				exponents.append(exp)
				if is_complex:
					generated_roots.append((int(root) ** exp))
				else:
					generated_roots.append((int(root) ** exp) % P)
		print('- All ' + str(M) + '-th roots of unity generated by the ' + str(exponents) + '-exponented primitive ' + str(M) + '-th root of unity ' + str(root) + ' : ' + str(generated_roots))
	print()




class Encoder:
	"""Basic BFV encoder to encode complex vectors into polynomials."""
	
	def __init__(self, N=None, P=None, scale=None, is_homogeneous_rotation=True):
		"""Initialization of the encoder for M a power of 2. 
		
		xi, which is an M-th root of unity will, be used as a basis for our computations.
		"""
		self.M = N*2
		self.N = N
		if P is not None:
			self.P = P
			self.N_inverse = mod_inv(N, P)
			self.is_complex=False
		else:
			self.N_inverse = 1/N
			self.is_complex=True
		self.scale = scale
		self.is_homogeneous_rotation = is_homogeneous_rotation
		self.create_sigma_R_basis()
		
	def vandermonde(self) -> np.array:
		"""Computes the Vandermonde matrix from a m-th root of unity."""
		matrix = []
		# We will generate each row of the matrix

		cyclic_x_values_ordered = []
		power_array = []
		power_matrix = []
		if self.is_complex:
			primitive_root_of_unity = np.exp(2 * np.pi * 1j / self.M)
		else:
			mult_group_generator = primitive_root(self.P)
			primitive_root_of_unity = pow(mult_group_generator, (self.P - 1) // self.M, self.P)

		for h in range(int(self.N/2)):
			new_root = primitive_root_of_unity ** JFunction(h, self.N, True)
			if not self.is_complex:
				new_root = centered_mod(new_root, self.P)
			if new_root in cyclic_x_values_ordered:
				 print("[Error] 1: root " + str(new_root) + " is already in " + str(cyclic_x_values_ordered))
				 sys.exit(0)
			cyclic_x_values_ordered.append(new_root)
			power_array.append(JFunction(h, self.N, True))
		for h in reversed(range(int(self.N/2))):
			new_root = primitive_root_of_unity ** JFunction(h, self.N, False)
			if not self.is_complex:
				new_root = centered_mod(new_root, self.P) 
			if new_root in cyclic_x_values_ordered:
				 print("[Error] 2: root " + str(new_root) + " is already in " + str(cyclic_x_values_ordered))
				 print("N: " + str(self.N))
				 sys.exit(0)
			cyclic_x_values_ordered.append(new_root)
			power_array.append(JFunction(h, self.N, False))

		for x in cyclic_x_values_ordered:
			# For each row we select a different root
			row = []
			# Then we store its powers
			for degree in range(self.N):
				insert = 1
				for count in range(degree):
					insert = insert * x if self.is_complex else centered_mod(insert * x, self.P)
				row.append(insert)
			matrix.append(row)
		return matrix

	def create_sigma_R_basis(self):
		"""Creates the basis (sigma(1), sigma(X), ..., sigma(X** N-1))."""
		self.sigma_R_basis = np.array(np.array(self.vandermonde()).T, dtype = object)

		# Convert each element to Python's built-in int type, otherwise the computation will overflow!
		#self.sigma_R_basis = MatrixDataSizeToInfinite(self.sigma_R_basis)
		self.sigma_R_basis_counter = self.sigma_R_basis.T

		if self.is_homogeneous_rotation:
			W_modified = copy.deepcopy(self.sigma_R_basis)
			WT_modified = copy.deepcopy(self.sigma_R_basis_counter)

			for i in range(self.N):
				for j in range(self.N // 4):
					val1 = W_modified[i][j]
					val2 = W_modified[i][self.N//2 - 1 - j]
					W_modified[i][j] = val2
					W_modified[i][self.N//2 - 1 - j] = val1

			for i in range(self.N//4):
					val1 = copy.deepcopy(WT_modified[self.N//2 + i])
					val2 = copy.deepcopy(WT_modified[self.N - 1 - i])
					WT_modified[self.N//2 + i] = val2
					WT_modified[self.N - 1 - i] = val1

			self.sigma_R_basis = W_modified
			self.sigma_R_basis_counter = WT_modified
		print()
		print("<W Matrix>")
		print(self.sigma_R_basis)
		print()
		print("<W^* Matrix>")
		print(self.sigma_R_basis_counter)
		print()
		print("<(W^*)*(W)>")
		if self.is_complex:
			print(np.matmul(self.sigma_R_basis_counter, self.sigma_R_basis))
		else:
			print(centered_mod_arr(np.matmul(self.sigma_R_basis_counter, self.sigma_R_basis), self.P))
		print()
					

	def encode(self, input_vector: np.array) -> Polynomial:
		#anti_diagonal_matrix = MatrixDataSizeToInfinite(np.eye(self.N)[::-1])
		anti_diagonal_matrix = np.eye(self.N, dtype=object)[::-1]
		if self.is_complex:
			for i in range(self.N // 2):
				lista = [input_vector[i].conjugate()]
				input_vector = np.append(input_vector, lista)
		basis_coordinates = np.matmul(self.N_inverse *  self.sigma_R_basis, anti_diagonal_matrix).dot(input_vector)
		if self.is_complex:
			basis_coordinates = np.asarray(basis_coordinates, dtype=object)
			basis_coordinates = np.array([c.real if isinstance(c, complex) else c for c in basis_coordinates], dtype=object)
			#basis_coordinates = np.real(basis_coordinates)
		else:
			basis_coordinates = centered_mod_arr(basis_coordinates, self.P)   # <-- centered mod P
		p = Polynomial(np.array(basis_coordinates, dtype=object))
		return p
		#scaled_p = p * self.scale
		#return scaled_p
	
	def decode(self, p: Polynomial) -> np.array:
		#rescaled_p = round_polynomial((p / self.scale))
		#coef = rescaled_p.coef
		#coef = p
		#coef_preserved_zeros = np.pad(coef, (0, self.N - coef.size), 'constant')
		coef = np.array(p.coef, dtype=object).reshape(-1)
		#coef_preserved_zeros = VectorDataSizeToInfinite(coef_preserved_zeros)
		z = np.matmul(self.sigma_R_basis_counter, coef)
		if not self.is_complex:
			z = centered_mod_arr(z, self.P)				   # <-- centered mod P
		return z



def JFunction(h, N, is_plus: bool):
	if is_plus: 
		return (5 ** h) % (2*N)
	else:
		return (-(5 ** h)) % (2*N)



'''
def VectorDataSizeToInfinite(vector):
	# Convert the matrix to dtype=object to allow arbitrary integer sizes
	vector_object = np.array(vector, dtype=object)
	print(vector)
	print(type(vector)) 
	print(type(vector.dtype)) 
	print(type(vector_object)) 
	print(type(vector_object.dtype)) 
	# Convert each element to Python's built-in int type
	for i in range(vector_object.shape[0]):
		vector_object[i] = int(vector_object[i])
	return vector_object
'''

def MatrixDataSizeToInfinite(matrix):
	# Convert the matrix to dtype=object to allow arbitrary integer sizes
	matrix_object = np.array(matrix, dtype=object)

	# Convert each element to Python's built-in int type
	for i in range(matrix_object.shape[0]):
		for j in range(matrix_object.shape[1]):
			matrix_object[i, j] = int(matrix_object[i, j])
	return matrix_object	



def round_polynomial(poly):
	degree = 0
	new_coef = []
	coef = poly.coef
	for co in coef:
		new_coef.append(round(co))
	poly.coef = new_coef
	return poly



def poly_mod(poly, N, Q=None):
	res = [0] * N
	for deg, co in enumerate(poly.coef):
		j = deg % N
		sign = -1 if ((deg // N) & 1) else 1   # x^n = -1 이므로 블록마다 부호 교대
		res[j] += sign * co
	if Q is not None:
		res = [centered_mod(r, Q) for r in res]	# 항상 centered residue
	return Polynomial(np.array(res, dtype=object))
	return poly



def poly_add(p1, p2, N, Q=None):
	return (poly_mod(p1 + p2, N, Q)) 



def poly_sub(p1, p2, N, Q=None):
	r = p1 - p2
	return (poly_mod(p1 - p2, N, Q)) 



def poly_mult(p1, p2, N, Q=None, is_ntt=True):
	if is_ntt and Q is not None:
		return poly_mult_ntt(p1, p2, N, Q)
	else:
		return poly_mult_naive(p1, p2, N, Q)



def get_ntt_root(M, Q):
	"""
	Finds a primitive M-th root of unity modulo Q.
	Requires Q = 1 mod M.
	"""
	if (Q - 1) % M != 0:
		return None
	
	# Q is prime, so primitive roots exist.
	# We try random candidates. g^( (Q-1)/M ) is a candidate.
	k = (Q - 1) // M
	for cand in range(2, Q):
		root = pow(cand, k, Q)
		# Check if root is primitive M-th root:
		# It must be that root^(M/2) != 1 (actually -1 mod Q)
		if pow(root, M // 2, Q) != 1:
			return root
	return None

def bit_reverse_copy(a, n):
	"""Permutes array 'a' of length 'n' by bit-reversing indices."""
	b = [0] * n
	bits = n.bit_length() - 1
	for i in range(n):
		# Bit reverse i
		rev = 0
		temp = i
		for _ in range(bits):
			rev = (rev << 1) | (temp & 1)
			temp >>= 1
		b[rev] = a[i]
	return b



def ntt_core(coeffs, N, Q, root, inverse=False):
	"""
	Iterative Cooley-Tukey NTT. 
	N must be a power of 2.
	"""
	n = len(coeffs)
	a = bit_reverse_copy(coeffs, n)
	
	# Precompute powers of root (twiddle factors)
	# If inverse, use modular inverse of root
	w_base = root
	if inverse:
		w_base = mod_inv(root, Q)
	
	m = 1
	while m < n:
		w_m = pow(w_base, n // (2 * m), Q)
		step = 2 * m
		for k in range(0, n, step):
			w = 1
			for j in range(m):
				t = (w * a[k + j + m]) % Q
				u = a[k + j]
				a[k + j] = (u + t) % Q
				a[k + j + m] = (u - t) % Q
				w = (w * w_m) % Q
		m = step
		
	if inverse:
		n_inv = mod_inv(n, Q)
		a = [(x * n_inv) % Q for x in a]
		
	return a



def poly_mult_ntt(p1, p2, N, Q):
	"""
	NTT-based multiplication.
	1. Zero-pad to 2*N (next power of 2).
	2. Forward NTT.
	3. Pointwise multiply.
	4. Inverse NTT.
	"""
	# 1. Setup size M (must be power of 2 >= 2N-1 for linear convolution)
	M = 1
	target_len = len(p1.coef) + len(p2.coef) - 1
	while M < target_len:
		M *= 2
		
	# 2. Check NTT suitability
	if not is_prime(Q) or (Q - 1) % M != 0:
		# Fallback to naive if Q is not an NTT-friendly prime
		return poly_mult_naive(p1, p2, N, Q)
		
	root = get_ntt_root(M, Q)
	if root is None:
		return poly_mult_naive(p1, p2, N, Q)

	# 3. Prepare coefficients (pad with zeros)
	c1 = [int(c) for c in p1.coef] + [0] * (M - len(p1.coef))
	c2 = [int(c) for c in p2.coef] + [0] * (M - len(p2.coef))
	
	# 4. NTT
	ntt_c1 = ntt_core(c1, M, Q, root, inverse=False)
	ntt_c2 = ntt_core(c2, M, Q, root, inverse=False)
	
	# 5. Pointwise Mult
	ntt_res = [(a * b) % Q for a, b in zip(ntt_c1, ntt_c2)]
	
	# 6. Inverse NTT
	res_coeffs = ntt_core(ntt_res, M, Q, root, inverse=True)
	
	# 7. Reduce X^N + 1 handled by poly_mod (truncating to N happens there)
	# We pass the full linear convolution result to poly_mod
	return poly_mod(Polynomial(np.array(res_coeffs, dtype=object)), N, Q)



def poly_mult_naive(p1, p2, N, Q=None):
	c1 = p1.coef
	c2 = p2.coef 
	conv = [0] * (len(c1) + len(c2) - 1)
	for i, a in enumerate(c1):
		if a == 0: continue
		for j, b in enumerate(c2):
			if b == 0: continue
			conv[i+j] += a*b
	return poly_mod(Polynomial(np.array(conv, dtype=object)), N, Q)



def poly_zero(N: int) -> Polynomial:
	return Polynomial(np.array([0] * N, dtype=object))



def poly_neg(p: Polynomial, N: int, Q: int) -> Polynomial:
	return poly_mod(Polynomial(np.array([-int(c) for c in p.coef], dtype=object)), N, Q)



def poly_scalar_mult(p: Polynomial, k: int, N: int, Q: int) -> Polynomial:
	kk = int(k)
	tmp = [int(c) * kk for c in p.coef]
	return poly_mod(Polynomial(np.array(tmp, dtype=object)), N, Q)



def poly_div_round(p: Polynomial, d: int, N: int, Q: int) -> Polynomial:
	# work on centered coefficients, then round-divide
	p = poly_mod(p, N, Q)
	coeffs = [round_div(int(c), int(d)) for c in p.coef]
	print("Coeff:", coeffs)
	return poly_mod(Polynomial(np.array(coeffs, dtype=object)), N, Q)



def sample_uniform_poly(n, q) -> Polynomial:
	return Polynomial(np.array([random.randrange(q) for _ in range(n)], dtype=object))



def sample_small_poly(n, bound = 1) -> Polynomial:
	return Polynomial(np.array([random.randint(-bound, bound) for _ in range(n)], dtype=object))



def digits_len(Q: int, base: int) -> int:
	"""
	Calculates L such that base^L > 2*Q.
	The factor of 2 is REQUIRED for Balanced Gadget Decomposition.
	"""
	Q = int(Q); base = int(base)
	L, cur = 0, 1
	# Check against 2 * Q, not just Q
	while cur < 2 * Q:
		cur *= base
		L += 1
	return max(L, 1)


def poly_rotate(poly, q, n, rotation_offset):
	degree = 0
	coef = poly.coef
	new_coef = [0] * n
	power = JFunction(rotation_offset, n, True)
	for i in range(len(coef)):
		if (((i * power) // n) % 2 == 1):
			new_coef[(i * power) % n] += -coef[i]
		else:
			new_coef[(i * power) % n] += coef[i]
	new_poly = Polynomial(np.array(list(new_coef), dtype=object))
	return new_poly


def poly_conjugate(poly, N, Q):
	"""
	Applies the automorphism X -> X^{-1} to the polynomial.
	Mathematically equivalent to X -> X^(2N-1).
	"""
	new_coef = [0] * N
	coef = poly.coef
	
	# The automorphism index for conjugation is 2N - 1
	power = 2 * N - 1
	
	for i in range(len(coef)):
		# Calculate new position: i * power
		# We need to handle the reduction modulo X^N + 1
		raw_deg = i * power
		new_deg = raw_deg % N
		
		# If (raw_deg // N) is odd, the term flips sign (because X^N = -1)
		if (raw_deg // N) % 2 == 1:
			new_coef[new_deg] += -coef[i]
		else:
			new_coef[new_deg] += coef[i]
			
	# Rebuild polynomial and ensure it is modulo Q
	return poly_mod(Polynomial(np.array(new_coef, dtype=object)), N, Q)



def gadget_decompose_poly(a: Polynomial, N: int, Q: int, base: int, L: int):
	"""
	Coefficient-wise gadget decomposition of a (mod Q):
		a ≡ Σ_{i=0..L-1} d_i * base^i   (mod Q)
	Returns list of polynomials d_i (dtype=object), using *balanced digits* to reduce noise.
	"""
	a = poly_mod(a, N, Q)
	coeffs = [int(c) for c in a.coef]

	digits = [[0] * N for _ in range(L)]
	halfB = base // 2

	# Coefficients beyond len(coeffs) are implicitly 0, so their digits remain 0.
	for j in range(min(N, len(coeffs))):
		x = coeffs[j]
		for i in range(L):
			d = x % base
			# balanced digit in roughly [-base/2, base/2]
			if d > halfB:
				d -= base
			digits[i][j] = int(d)
			x = (x - d) // base

	return [Polynomial(np.array(digits[i], dtype=object)) for i in range(L)]

	
	
	
def key_generate(n) -> Polynomial:
	"""Binary secret key s in {-1, 0, 1}^n """
	return Polynomial(np.array([random.randint(-1, 1) for _ in range(n)], dtype=object))
	
	
def fhe_keyswitch(a: Polynomial, b: Polynomial, ksk_bundle):
	"""
	Input ciphertext (a,b) under s_old.
	Output ciphertext (a2,b2) under s_new, using precomputed ksk_bundle from fhe_keyswitch_keygen().

	Uses:
		b2 = b + Σ d_i * B_i
		a2 = - Σ d_i * A_i
	where {d_i} is gadget decomposition of 'a'.
	"""
	base = int(ksk_bundle["base"])
	L	= int(ksk_bundle["L"])
	N	= int(ksk_bundle["N"])
	Q	= int(ksk_bundle["Q"])
	ksk  = ksk_bundle["ksk"]

	digits = gadget_decompose_poly(a, N, Q, base, L)

	accA = poly_zero(N)
	accB = poly_zero(N)

	for i in range(L):
		di = digits[i]
		Ai, Bi = ksk[i]
		accA = poly_add(accA, poly_mult(di, Ai, N, Q), N, Q)
		accB = poly_add(accB, poly_mult(di, Bi, N, Q), N, Q)

	a2 = poly_neg(accA, N, Q)
	b2 = poly_add(b, accB, N, Q)
	return (a2, b2)


# ---------- key-switch key precomputation ----------
def fhe_keyswitch_keygen(s_old: Polynomial, s_new: Polynomial, N, Q, base: int = 1 << 10, bound: int = 1, bgv_scale: int = None):
	"""
	Build KSK for switching from old secret key s_old to new secret key s_new.

	For each i, KSK[i] = (A_i, B_i) such that:
		B_i + A_i * s_new  ≈  (base^i) * s_old   (mod Q)
	"""
	L = digits_len(Q, base)
	ksk = []
	base_pow = 1

	for i in range(L):
		# message to embed (NO scale here): base^i * s_old  in R_Q
		msg_i = poly_scalar_mult(s_old, base_pow, N, Q)

		A_i = sample_uniform_poly(N, Q)			  # already dtype=object
		E_i = sample_small_poly(N, bound=bound)

		if bgv_scale is not None: # in the case of BGV
			E_i = poly_scalar_mult(E_i, int(bgv_scale), N, Q)

		tmp = poly_mult(A_i, s_new, N, Q)
		tmp = poly_add(tmp, E_i, N, Q)
		B_i = poly_add(tmp, msg_i, N, Q)

		ksk.append((A_i, B_i))
		base_pow *= base

	return {"base": base, "L": L, "N": N, "Q": Q, "ksk": ksk}
	
