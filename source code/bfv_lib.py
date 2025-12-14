import numpy as np
from numpy.polynomial import Polynomial
from fastcore.foundation import patch_to
from os import sys
import random
import math
from sympy import nextprime, mod_inverse, primitive_root
import copy



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

def EvaluateRandomCyclotomicPolynomial(x):
	return (x ** N) + 1

def IsPowerOf2(N):
	n_checker = N
	while n_checker > 1:
		if n_checker % 2 == 1:
			return False
		n_checker /= 2 
	return True


def Egcd(a, b):
	if a == 0:
		return (b, 0, 1)
	else:
		g, y, x = Egcd(b % a, a)
		return (g, x - (b // a) * y, y)

def ModInv(a, m):
	g, x, y = Egcd(a, m)
	if g != 1:
		#raise Exception('modular inverse does not exist')
		return -1
	else:
		return x % m


np.set_printoptions(suppress=True)


def JFunction(h, N, is_plus: bool):
	if is_plus: 
		return (5 ** h) % (2*N)
	else:
		return (-(5 ** h)) % (2*N)

def isprime(n):
	'''check if integer n is a prime'''
	# make sure n is a positive integer
	n = abs(int(n))
	# 0 and 1 are not primes
	if n < 2:
		return False
	# 2 is the only even prime number
	if n == 2: 
		return True	
	# all other even numbers are not primes
	if not n & 1: 
		return False
	# range starts with 3 and only needs to go up the squareroot of n
	# for all odd numbers
	for x in range(3, int(n**0.5)+1, 2):
		if n % x == 0:
			return False
	return True

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

class BFVEncoder:
	"""Basic BFV encoder to encode complex vectors into polynomials."""
	
	def __init__(self, P, Q, N, is_homogeneous_rotation):
		"""Initialization of the encoder for M a power of 2. 
		
		xi, which is an M-th root of unity will, be used as a basis for our computations.
		"""
		self.M = N*2
		self.P = P
		self.Q = Q
		self.scale = math.floor(Q / P)
		self.N = N
		self.N_inverse = ModInv(N, P)
		self.is_homogeneous_rotation = is_homogeneous_rotation
		self.create_sigma_R_basis()
		
	def vandermonde(self) -> np.array:
		"""Computes the Vandermonde matrix from a m-th root of unity."""
		matrix = []
		# We will generate each row of the matrix

		cyclic_x_values_ordered = []
		power_array = []
		power_matrix = []
		mult_group_generator = primitive_root(self.P)
		primitive_root_of_unity = pow(mult_group_generator, (self.P - 1) // self.M, self.P)
		for h in range(int(self.N/2)):
			new_root = centered_mod(primitive_root_of_unity ** JFunction(h, self.N, True), self.P)
			if new_root in cyclic_x_values_ordered:
				 print("Error 1: root " + str(new_root) + " is already in " + str(cyclic_x_values_ordered))
				 sys.exit(0)
			cyclic_x_values_ordered.append(new_root)
			power_array.append(JFunction(h, self.N, True))
		for h in reversed(range(int(self.N/2))):
			new_root = centered_mod((primitive_root_of_unity ** JFunction(h, self.N, False)), self.P) 
			if new_root in cyclic_x_values_ordered:
				 print("Error 2: root " + str(new_root) + " is already in " + str(cyclic_x_values_ordered))
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
					insert = centered_mod(insert * x, self.P)
				row.append(insert)
			matrix.append(row)
		return matrix

	def create_sigma_R_basis(self):
		"""Creates the basis (sigma(1), sigma(X), ..., sigma(X** N-1))."""
		self.sigma_R_basis = np.array(np.array(self.vandermonde()).T, dtype = object)

		# Convert each element to Python's built-in int type, otherwise the computation will overflow!
		self.sigma_R_basis = MatrixDataSizeToInfinite(self.sigma_R_basis)
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
		print(centered_mod_arr(np.matmul(self.sigma_R_basis_counter, self.sigma_R_basis), self.P))
		print()
					

	def encode(self, input_vector: np.array) -> Polynomial:
		anti_diagonal_matrix = MatrixDataSizeToInfinite(np.eye(self.N)[::-1])
		raw = np.matmul(self.N_inverse * self.sigma_R_basis, anti_diagonal_matrix).dot(input_vector)
		basis_coordinates = centered_mod_arr(raw, self.P)   # <-- centered mod P
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
		raw = np.matmul(self.sigma_R_basis_counter, coef)
		z = centered_mod_arr(raw, self.P)				   # <-- centered mod P
		return z

def round_polynomial(poly):
	degree = 0
	new_coef = []
	coef = poly.coef
	for co in coef:
		new_coef.append(round(co))
	poly.coef = new_coef
	return poly


def poly_mod(poly, n, q):
	res = [0] * n
	for deg, co in enumerate(poly.coef):
		j = deg % n
		sign = -1 if ((deg // n) & 1) else 1   # x^n = -1 이므로 블록마다 부호 교대
		res[j] += sign * co

	res = [centered_mod(r, q) for r in res]	# 항상 centered residue
	return Polynomial(np.array(res, dtype=object))
	return poly

def poly_add(p1, p2, n, q):
	return (poly_mod(p1 + p2, n, q)) 

def poly_sub(p1, p2, N, Q):
	r = p1 - p2
	return (poly_mod(p1 - p2, N, Q)) 

def poly_mult(p1, p2, n, q):
	#p1 = as_int_obj_poly(p1); p2 = as_int_obj_poly(p2)
	#c1 = list(map(p1.coef)); c2 = list(map(p2.coef))
	c1 = p1.coef
	c2 = p2.coef 
	#if not c1 or not c2:
	#	print("Error: poly_mult input invalid")
	#	print(p1, p2, c1, c2)
	#	sys.exit(0)
		#return Polynomial(np.array([0]*n, dtype=object))
	conv = [0] * (len(c1) + len(c2) - 1)
	for i, a in enumerate(c1):
		if a == 0: continue
		for j, b in enumerate(c2):
			if b == 0: continue
			conv[i+j] += a*b
	return poly_mod(Polynomial(np.array(conv, dtype=object)), n, q)

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
	"""Small integer L such that base^L >= Q."""
	Q = int(Q); base = int(base)
	L, cur = 0, 1
	while cur < Q:
		cur *= base
		L += 1
	return max(L, 1)

def gadget_decompose_poly(a: Polynomial, N: int, Q: int, base: int, L: int):
	"""
	Coefficient-wise gadget decomposition of a (mod Q):
		a ≡ Σ_{i=0..L-1} d_i * base^i   (mod Q)
	Returns list of polynomials d_i (dtype=object), using *balanced digits* to reduce noise.
	"""
	coeffs = [int(c) % Q for c in a.coef]			  # 0..Q-1 representatives

	digits = [[0] * N for _ in range(L)]
	halfB = base // 2

	for j in range(N):
		x = coeffs[j]
		for i in range(L):
			d = x % base
			# balanced digit in roughly [-base/2, base/2]
			if d > halfB:
				d -= base
			digits[i][j] = int(d)
			x = (x - d) // base

	return [Polynomial(np.array(digits[i], dtype=object)) for i in range(L)]


# ---------- key-switch key precomputation ----------
def fhe_keyswitch_keygen(s_old: Polynomial, s_new: Polynomial, N, Q, base: int = 1 << 10, bound: int = 1):
	"""
	Build KSK for switching from old secret key s_old to new secret key s_new.

	For each i, KSK[i] = (A_i, B_i) such that:
		B_i - A_i * s_new  ≈  (base^i) * s_old   (mod Q)
	(same "minus" convention: b - a*s)
	"""
	L = digits_len(Q, base)
	ksk = []
	base_pow = 1

	for i in range(L):
		# message to embed (NO scale here): base^i * s_old  in R_Q
		msg_i = poly_scalar_mult(s_old, base_pow, N, Q)

		A_i = sample_uniform_poly(N, Q)			  # already dtype=object
		E_i = sample_small_poly(N, bound=bound)

		tmp = poly_mult(A_i, s_new, N, Q)
		tmp = poly_add(tmp, E_i, N, Q)
		B_i = poly_add(tmp, msg_i, N, Q)

		ksk.append((A_i, B_i))
		base_pow *= base

	return {"base": base, "L": L, "N": N, "Q": Q, "ksk": ksk}

def fhe_relin_keygen(s: Polynomial, N: int, Q: int, base: int = 1 << 10, bound: int = 1):
	"""
	Build relin key for ct×ct multiplication:
	  RLK = KeySwitchKey(s_old = s^2, s_new = s)

	It returns the same bundle format as fhe_keyswitch_keygen().
	"""
	s2 = poly_mult(s, s, N, Q)
	return fhe_keyswitch_keygen(s2, s, N, Q, base=base, bound=bound)


def key_generate(n) -> Polynomial:
	"""Binary secret key s in {-1, 0, 1}^n """
	global _BFV_SK
	_BFV_SK = Polynomial(np.array([random.randint(-1, 1) for _ in range(n)], dtype=object))
	return _BFV_SK


def fhe_encrypt(m: Polynomial, N, P, Q, s, scale):
	"""
	Returns ciphertext (a, b):
	  a <- uniform in R_q
	  e <- small
	  b = a*s + e + scale*m   (mod x^n+1, q)
	"""
	m_red = poly_mod(m, N, P)  # coeffs in [0,p)
	m_scaled = Polynomial(np.array([c * scale for c in m_red.coef], dtype=object))

	a = sample_uniform_poly(N, Q)
	e = sample_small_poly(N, bound=1)

	b = poly_add(poly_add(poly_mult(a, s, N, Q), e, N, Q),
				 m_scaled, N, Q)

	# store centered representatives (optional but consistent with your “always centered” style)
	a = poly_mod(a, N, Q)
	b = poly_mod(b, N, Q)
	return (a, b)

def fhe_decrypt(a: Polynomial, b: Polynomial, N, P, Q, s, scale) -> Polynomial:
	"""
	Decrypt:
	  v = b - a*s  (mod x^n+1, q)  ~= scale*m + e
	  m_hat = round(v/scale)
	  (optionally centered mod p if p was stored by fhe_encrypt)
	"""

	v = poly_sub(b, poly_mult(a, s, N, Q), N, Q)
	m_hat = [round_div(c, scale) for c in v.coef]
	temp = poly_mult(a, s, N, Q)
	m_hat = [centered_mod(c, P) for c in m_hat]

	return Polynomial(np.array(m_hat, dtype=object))

def fhe_add_cipher_cipher(a1, b1, a2, b2, N, Q):
	a3 = poly_add(a1, a2, N, Q)
	b3 = poly_add(b1, b2, N, Q)
	return (a3, b3)

def fhe_add_cipher_plain(a1, b1, p2, N, Q, scale):
	b3 = poly_add(b1, poly_scalar_mult(p2, scale, N, Q), N, Q)
	return (a1, b3) 

def fhe_mult_cipher_plain(a1, b1, p2, q, n):
	a3 = poly_mult(a1, p2)
	b3 = poly_mult(b1, p2) 
	return (a3, b3)

def fhe_rotate(a1, b1, N, Q, rotation_offset, galois_rotation_key):
    # 1) Automorphism on ciphertext components
    a_rot = poly_mod(poly_rotate(a1, Q, N, rotation_offset), N, Q)
    b_rot = poly_mod(poly_rotate(b1, Q, N, rotation_offset), N, Q)

    # 2) Key-switch from s_rot -> s (cache by (offset,N,Q,s)) by gadget decomposition
    a_out, b_out = fhe_keyswitch(a_rot, b_rot, galois_rotation_key)
    return (a_out, b_out)

	

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

def fhe_keyswitch(a: Polynomial, b: Polynomial, ksk_bundle):
	"""
	Input ciphertext (a,b) under s_old.
	Output ciphertext (a2,b2) under s_new, using precomputed ksk_bundle from fhe_keyswitch_keygen().

	Uses:
		b2 = b - Σ d_i * B_i
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
	b2 = poly_sub(b, accB, N, Q)
	return (a2, b2)

def poly_rescale_P_over_Q(raw_poly: Polynomial, P: int, Q: int, N: int) -> Polynomial:
    """
    BFV rescale: round((P/Q) * coeff) coefficient-wise, then reduce mod Q.
    Output is in R_Q.
    """
    raw = [int(x) for x in np.array(raw_poly.coef, dtype=object).reshape(-1)]
    scaled = [round_div(P * c, Q) for c in raw]     # big-int rounding
    return poly_mod(Polynomial(np.array(scaled, dtype=object)), N, Q)

def fhe_mult_cipher_cipher(a1, b1, a2, b2, N, P, Q, scale, relin_key):
    """
    BFV ct×ct with gadget decomposition relinearization.
    IMPORTANT: uses BFV rescale round(P/Q * ·), NOT division by 'scale'.
    """
    # Use the plaintext modulus from the relin key bundle if you store it elsewhere.
    # Here we rely on the module global P being passed through your encrypt/decrypt code path.
    # We'll fetch it from relin_key if you add it, otherwise require a global.

    # Ensure inputs are centered in R_Q
    a1 = poly_mod(a1, N, Q); b1 = poly_mod(b1, N, Q)
    a2 = poly_mod(a2, N, Q); b2 = poly_mod(b2, N, Q)

    # ---- Raw products in big integers (no mod Q yet) ----
    d0_raw = poly_mult(b1, b2, N, Q*scale)  # b1*b2
    d1_raw = poly_add(
        poly_mult(b1, a2, N, Q*scale),
        poly_mult(a1, b2, N, Q*scale),
        N, Q*scale
    )  # b1*a2 + a1*b2
    d2_raw = poly_mult(a1, a2, N, Q*scale)  # a1*a2

    # ---- BFV rescale: round(P/Q * ·) and reduce mod Q ----
    d0 = poly_rescale_P_over_Q(d0_raw, P, Q, N)
    d1 = poly_rescale_P_over_Q(d1_raw, P, Q, N)
    d2 = poly_rescale_P_over_Q(d2_raw, P, Q, N)

    # ---- Relinearize the s^2 term using your gadget-decomp keyswitch ----
    # Under s^2, (a_old=-d2, b_old=0) has phase = d2*s^2
    a_old = poly_mod(poly_neg(d2, N, Q), N, Q)
    b_old = poly_zero(N)
    a_sw, b_sw = fhe_keyswitch(a_old, b_old, relin_key)

    # Combine into 2-term ciphertext under s
    a3 = poly_add(d1, a_sw, N, Q)
    b3 = poly_add(d0, b_sw, N, Q)
    return (a3, b3)


def fhe_bootstrapping(a1, b1, N, P, Q):
	return 0
