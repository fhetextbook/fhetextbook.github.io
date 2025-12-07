import numpy as np
from numpy.polynomial import Polynomial
from fastcore.foundation import patch_to
from os import sys
import random
import math
from sympy import nextprime, mod_inverse, primitive_root
import copy

# First we set the parameters

#is_random_param = True
is_random_param = False
is_homogeneous_rotation = True

#############################
# Test Option 1: Targeted Test Params


# Testcase 1
#N = 4
#P = 17
#bfv_vector = [10, 3, 5, 13]
#bfv_vector2 = [2, 4, 3, 6]
#rotation_offset = 3


# Testcase 2
N = 4
P = 17
bfv_vector = [10, 3, 5, 13]
bfv_vector2 = [2, 4, 3, 6]
rotation_offset = 3


Q = 2**8
scale = math.floor(Q / P)
print('scale : ' + str(scale))

def EvaluateTargetedCyclotomicPolynomial(x):
	return (x ** N) + 1

#############################
# Test Option 2: Random Test Params

test_count = 1000

N_LOWER = 8
N_UPPER = 8

P_LOWER = 35
P_UPPER = 10000

scale_LOWER = 1
scale_UPPER = 64

def EvaluateRandomCyclotomicPolynomial(x):
	return (x ** N) + 1

##############################################



N_inverse = -1
primitive_root_of_unity = -1

test_index = 0
M = N*2

def IsPowerOf2(N : int):
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


def GenerateParam():
	global N
	global P
	global bfv_vector
	global gfv_vector2
	global bfv_vector
	global bfv_vector2
	global M
	global N_inverse
	global primitive_root_of_unity
	#global scale_inverse
	global scale
	global Q

	if is_random_param:
		# Decide the exponent N
		while True:
			power = random.randint(N_LOWER, N_UPPER)
			if IsPowerOf2(power):
				 N = power
				 P = random.randint(P_LOWER, P_UPPER)
				 #scale = random.randint(scale_LOWER, scale_UPPER)
				 M = N * 2
				 break
			rotation_offset = random.randint(0, N // 2 - 1)
	else:
		if len(bfv_vector) != N or len(bfv_vector2) != N:
		  print("Error: The dimensions of the input vectors are not " + str(N) + ", but " + str(len(bfv_vector)) + ", " + str(len(bfv_vector2)))
		  sys.exit()
	if M <= 0:
		print("Error: M has to be specified for a targeted test")
		sys.exit()

	# Fermit's Little theorem: if A and p are co-primes, then A^{p-1} ≡ 1 mod p 
	# ==> If p is a prime, then any number A in [0, p - 1] has co-prime with p, thus any number A in [0, p - 1] has a multiplicative inverse, which is A^{p-2} mod p, and A * A^{p-2} = A^{p-1} ≡ 1 mod p

	# Find a modulus P that guarnatees the existence of a primitive M-th root of unity and the inverse of N
	# If P is a prime, then the inverse of N always exists
	# A primitive (M=2N)-th root of unity exists if M=2N divides P - 1 ???
	#  X is the primitive M-th root of unity if X^M ≡ 1 mod P AND Order(X) = M (where M = 2N)

	if is_random_param:
		P = nextprime(M)  # Find a prime P bigger than M(=2N)
		# We find P such that there exists the primitive M-th root of unity in modulo P. According to Fermat's Little Theorem (A^{P-1} ≡ 1 mod P), each element in P has an order that divides P-1 (i.e., the multiplicative modulo group size), and all elements in the group collectively cover all divisors of P-1. Among these elements, given Ord(A) = M for some element A, that's the primitive M-th root of unity modulo P (as the definition of primitive M-th root of unity: any M-th root of unity C, solution for X^M ≡ 1, such that Ord(C) = M). Therefore, if P-1 (the multiplicative modulo group size) is some multiple of M, then A^{P-1} = A^{M*k} = (A^k)^M ≡ 1 mod P, then some element B=A^k in this group will have Ord(B) = M. As a conclusion, if (P - 1) mod M == 0 (i.e., M divides the size of the multiplicative group of P), then some elements in this group will have an order M, which will be the primitive M-th root of unity. 
		while (P - 1) % M != 0: 
			P = nextprime(P)
	print("Final P: " + str(P))
	N_inverse = ModInv(N, P) # If P is a prime, the multiplicative inverse of N (< P) is guaranteed to exist (by Fermit's Little theorem)
	#Q = P*scale*scale
	#print(Q)
	#print(scale)
	root_values = []

	# [OPTION 1]: A brute-force way of finding the primitiev (M=2N)-th root of unity (i.e., the solution for X^N + 1 = 0 mod P)
	'''
	# After then, we find X ∈ {1,2,...P-1} whose order is 2N. That value is the primitive P-th root of unity.
	while len(root_values) == 0:

		# Find P such that X^{P-1} ≡ 1 mod P  (Fermat's little theorem) AND M=2N divides P-1, which ensures that the order of each X ∈ {1,2,...P-1} mod P is either 2N or some factor of 2N

		N_inverse = ModInv(N, P) # If P is a prime, the multiplicative inverse of N (< P) is guaranteed to exist (by Fermit's Little theorem)
		if not is_random_param and N_inverse <= 0:
			print("Error: N=" + str(N) + "has no inverse modulo P=" + str(P))
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
					root_values.append(i)   # Then, i is the primitive M-th root of unity
		#if not is_random_param: 
		if len(root_values) == 0 and not is_random_param: 
			print("Modulus P=" + str(P) + " has no primitive " + str(M) + "-th root of unity")
			sys.exit(0)

	primitive_root_of_unity = min(root_values)
	'''

#	'''
	### [OPTION 2]: An alternative efficient way of computing the primitive (M=2N)-th root of unity
	mult_group_generator = primitive_root(P)  # prim_root is a number such that Order(prim_root) = totient(p) mod p
	primitive_root_of_unity = pow(mult_group_generator, (P - 1) // M, P) # The (M=2N)-th primitive root of unity (i.e., solution for both X^N + 1 = 0 mod p AND X^{2N} = 1 mod P)
	print(primitive_root_of_unity)
	next_root = primitive_root_of_unity
	while True:
		root_values.append(next_root)
		next_root = next_root * (primitive_root_of_unity ** 2) % P
		if next_root == primitive_root_of_unity:
			break
#	'''

	print("==================== PARAMETERS ====================")
	print("N: " + str(N))
	print("Plaintext modulus: " + str(P))
	print("Ciphertext modulus: " + str(Q))
	print("Scale Factor: " + str(scale))
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
				generated_roots.append((int(root) ** exp) % P)
		print('- Roots of the ' + str(M) + '-th cyclotomic polynomial generated by the co-prime-' + str(coprimes) + '-exponented primitive ' + str(M) + '-th root of unity ' + str(root) + ' : ' + str(generated_roots))
	print()
	for root in root_values:
		generated_roots = []
		exponents = []
		for exp in range(M):
				exponents.append(exp)
				generated_roots.append((int(root) ** exp) % P)
		print('- All ' + str(M) + '-th roots of unity generated by the ' + str(exponents) + '-exponented primitive ' + str(M) + '-th root of unity ' + str(root) + ' : ' + str(generated_roots))
	print()

	if is_random_param:
		bfv_vector = []
		bfv_vector2 = []
		for i in range(N):
			bfv_vector.append(random.randint(0, P - 1))
			bfv_vector2.append(random.randint(0, P - 1))

np.set_printoptions(suppress=True)


#print("- Polynomial Ring: X^" + str(N) + " + 1")
#print()


def JFunction(h, is_plus: bool):
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


def VectorDataSizeToInfinite(vector):
	# Convert the matrix to dtype=object to allow arbitrary integer sizes
	vector_object = np.array(vector, dtype=object)

	# Convert each element to Python's built-in int type
	for i in range(vector_object.shape[0]):
		vector_object[i] = int(vector_object[i])
	return vector_object

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
	
	def __init__(self):
		"""Initialization of the encoder for M a power of 2. 
		
		xi, which is an M-th root of unity will, be used as a basis for our computations.
		"""
		self.M = M
		self.create_sigma_R_basis()
		

	@staticmethod
	def vandermonde(M: int) -> np.array:
		"""Computes the Vandermonde matrix from a m-th root of unity."""
		global P	 
		global N_inverse
		global N
		matrix = []
		# We will generate each row of the matrix

		cyclic_x_values_ordered = []
		power_array = []
		power_matrix = []
		for h in range(int(N/2)):
			new_root = (primitive_root_of_unity ** JFunction(h, True)) % P 
			if new_root in cyclic_x_values_ordered:
				 print("Error 1: root " + str(new_root) + " is already in " + str(cyclic_x_values_ordered))
				 sys.exit(0)
			cyclic_x_values_ordered.append(new_root)
			power_array.append(JFunction(h, True))
		for h in reversed(range(int(N/2))):
			new_root = (primitive_root_of_unity ** JFunction(h, False)) % P 
			if new_root in cyclic_x_values_ordered:
				 print("Error 2: root " + str(new_root) + " is already in " + str(cyclic_x_values_ordered))
				 print("N: " + str(N))
				 sys.exit(0)
			cyclic_x_values_ordered.append(new_root)
			power_array.append(JFunction(h, False))

		for x in cyclic_x_values_ordered:
			# For each row we select a different root
			row = []
			# Then we store its powers
			for degree in range(N):
				insert = 1
				for count in range(degree):
					insert = insert * x % P
				row.append(insert)
			matrix.append(row)
		#for x in power_array:
			# For each row we select a different root
			#row = []
			# Then we store its powers
			#for degree in range(N):
			#	row.append(x ** degree)
			#power_matrix.append(row)
		return matrix

	def create_sigma_R_basis(self):
		"""Creates the basis (sigma(1), sigma(X), ..., sigma(X** N-1))."""
		self.sigma_R_basis = np.array(np.array(self.vandermonde(self.M)).T, dtype = object)

		# Convert each element to Python's built-in int type, otherwise the computation will overflow!
		self.sigma_R_basis = MatrixDataSizeToInfinite(self.sigma_R_basis)
		self.sigma_R_basis_counter = self.sigma_R_basis.T

		if is_homogeneous_rotation:

			W_modified = copy.deepcopy(self.sigma_R_basis)
			WT_modified = copy.deepcopy(self.sigma_R_basis_counter)

			for i in range(N):
				for j in range(N // 4):
					val1 = W_modified[i][j]
					val2 = W_modified[i][N//2 - 1 - j]
					W_modified[i][j] = val2
					W_modified[i][N//2 - 1 - j] = val1

			for i in range(N//4):
					val1 = copy.deepcopy(WT_modified[N//2 + i])
					val2 = copy.deepcopy(WT_modified[N - 1 - i])
					WT_modified[N//2 + i] = val2
					WT_modified[N - 1 - i] = val1

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
		print(np.matmul(self.sigma_R_basis_counter, self.sigma_R_basis) % P)
		print()
					
@patch_to(BFVEncoder)
def __init__(self):
	self.M = M
	self.create_sigma_R_basis()
	self.scale = scale
	
@patch_to(BFVEncoder)
def encode(self, input_vector: np.array) -> Polynomial:
	anti_diagonal_matrix = MatrixDataSizeToInfinite(np.eye(N)[::-1])
	#print(self.sigma_R_basis)
	#print(anti_diagonal_matrix)
	#print(scaled_input_vector)
	basis_coordinates = np.matmul(N_inverse *  self.sigma_R_basis, anti_diagonal_matrix).dot(input_vector) % P
	p = Polynomial(basis_coordinates)
	#print("Inverse: " + str(N_inverse))
	#print("Input Vector: " + str(input_vector))
	#print("Poly: " + str(p.coef))
	scaled_p = p * self.scale
	#print("Scale: " + str(self.scale))
	#print("Check")
	#for element in self.sigma_R_basis:
	#	for element2 in element:
	#		print(f"Element: {element2}, Type: {type(element2)}")
	#for element in self.sigma_R_basis.T:
	#	for element2 in element:
	#		print(f"Element: {element2}, Type: {type(element2)}")
	#print(self.sigma_R_basis)
	#print(self.sigma_R_basis.T)
	#print(np.matmul(self.sigma_R_basis.T, self.sigma_R_basis) % P)
	#print(self.sigma_R_basis.T % P)
	#print(self.sigma_R_basis % P)
	#print()
	#print(self.sigma_R_basis.T[2])
	#print(self.sigma_R_basis.T[N - 3])
	#print()
	#print(self.sigma_R_basis.T[0] % P)
	#print(self.sigma_R_basis.T[N - 1] % P)
	#print()
	#arr1 = self.sigma_R_basis.T[0]
	#arr2 = self.sigma_R_basis.T[N - 1]
	#for element in arr1:
	#	print(f"Element: {element}, Type: {type(element)}")
	#print(np.dot(arr1, arr2))
	#total = 0
	#for i in range(N):
	#	total += arr1[i] * arr2[i]
	#print(total)
	#print()
	#print(self.sigma_R_basis.T[0].dot(self.sigma_R_basis.T[N - 1]) % P)
	#print((self.sigma_R_basis.T[0] % P).dot(self.sigma_R_basis.T[N - 1] % P) % P)
	#print()
	#print(np.matmul(self.sigma_R_basis.T, self.sigma_R_basis) % P)
	#print(np.matmul(self.sigma_R_basis.T % P, self.sigma_R_basis % P) % P)
	#print()
	#print(np.matmul(self.sigma_R_basis.T, self.sigma_R_basis))
	#print(anti_diagonal_matrix)
	#a = np.matmul(anti_diagonal_matrix, np.matmul(self.sigma_R_basis.T, self.sigma_R_basis))
	#print(a)
	#print(a % P)
	#print(np.matmul(self.sigma_R_basis.T, self.sigma_R_basis) % P)
	#print(np.matmul(anti_diagonal_matrix, np.matmul(self.sigma_R_basis.T, self.sigma_R_basis) % P))
	#print(np.matmul(anti_diagonal_matrix, np.matmul(self.sigma_R_basis.T, self.sigma_R_basis)) % P)
	#print()
	#print(np.matmul(np.eye(N)[::-1], np.matmul(self.sigma_R_basis.T, self.sigma_R_basis)).dot(scaled_input_vector))
	#print(np.matmul(np.eye(N)[::-1], np.matmul(self.sigma_R_basis.T, self.sigma_R_basis)).dot(scaled_input_vector) % P)
	#print((np.matmul(np.eye(N)[::-1], np.matmul(self.sigma_R_basis.T, self.sigma_R_basis))) % P)
	#print((N_inverse * np.matmul(np.eye(N)[::-1], np.matmul(self.sigma_R_basis.T, self.sigma_R_basis))) % P)
	#print()
	#print(scaled_input_vector)
	#print((np.matmul(np.eye(N)[::-1], np.matmul(self.sigma_R_basis.T, self.sigma_R_basis)) % P).dot(scaled_input_vector))
	#print(((np.matmul(np.eye(N)[::-1], np.matmul(self.sigma_R_basis.T, self.sigma_R_basis))).dot(scaled_input_vector)) % P)
	#print((np.matmul(np.eye(N)[::-1], np.matmul(self.sigma_R_basis.T, self.sigma_R_basis)) % P).dot(scaled_input_vector) % P)
	#print(N_inverse * (np.matmul(np.eye(N)[::-1], np.matmul(self.sigma_R_basis.T, self.sigma_R_basis))).dot(scaled_input_vector) % P)
	#print(np.matmul(self.sigma_R_basis.T, np.matmul(N_inverse *  self.sigma_R_basis, np.eye(N)[::-1])).dot(scaled_input_vector) % P)
	#print(np.matmul(self.sigma_R_basis.T, (np.matmul(N_inverse *  self.sigma_R_basis, np.eye(N)[::-1])).dot(scaled_input_vector)) % P)
	#print()

	return scaled_p

@patch_to(BFVEncoder)
def decode(self, p: Polynomial) -> np.array:
	rescaled_p = round_polynomial((p / scale))
	#print(self.sigma_R_basis)
	#print(rescaled_p.coef)
	coef = rescaled_p.coef
	coef_preserved_zeros = np.pad(coef, (0, N - len(coef)), 'constant')
	coef_preserved_zeros = VectorDataSizeToInfinite(coef_preserved_zeros)
	z = np.matmul(self.sigma_R_basis_counter, coef_preserved_zeros) % P
	#print("R Basis")
	#print(self.sigma_R_basis.T)
	#print("Coefficiet Mult")
	#print(coef_preserved_zeros)
	#print("Multiply")
	#print(z)
	#print()

	return z


def round_polynomial(poly):
	degree = 0
	new_coef = []
	coef = poly.coef
	for co in coef:
		new_coef.append(round(co))
	poly.coef = new_coef
	return poly

def reduce_polynomial(poly, n, q):
	degree = 0
	new_coef = []
	coef = poly.coef
	for co in coef:
		if degree >= n:
			new_coef[degree % n] += co * (-1 ** (degree / n))
			new_coef[degree % n] %= q 
		else:
			new_coef.append(co % q)
		degree += 1
	poly.coef = new_coef
	return poly

def rotate_homomorphic(poly, rotation_offset, n):
	degree = 0
	coef = poly.coef
	new_coef = [0] * n
	#new_coef = VectorDataSizeToInfinite(new_coef)
	power = JFunction(rotation_offset, True)
	#print("Power: " + str(power))
	for i in range(len(coef)):
		#print("Powered Coeff " + str(i) + " : " + str(i * power) + " --> " + str((i * power) % n)) 
		#print((i * power) % n)
		#print(new_coef)
		#print(len(coef))
		if (((i * power) // n) % 2 == 1):
			new_coef[(i * power) % n] += -coef[i]
		else:
			new_coef[(i * power) % n] += coef[i]
	poly.coef = new_coef
	return poly


def main():

	print("======================================")
	global M
	global N_inverse
	global primitive_root_of_unity
	global scale

	GenerateParam()
	encoder = BFVEncoder()

	print("<The BFV plaintext vector to encode>")
	z = np.array(bfv_vector)
	z2 = np.array(bfv_vector2)
	z3 = (z + z2) % P
	z4 = z * z2
	print('vector 1: ' + str(z))
	#p = encoder.encode(z)
	#decoded_z = encoder.decode(p)
	#print('decoded vector 1: ' + str(decoded_z))
	#sys.exit()
	print('vector 2: ' + str(z2))
	print('vector 1 + 2: ' + str(z3 % P))
	print('vector 1 * 2: ' + str(z4 % P))
	print()

	print("<The encoded BFV plaintext polynomial>")
	p = encoder.encode(z)
	p2 = encoder.encode(z2)
	p3 = reduce_polynomial(p + p2, N, Q)
	p3_rotated = reduce_polynomial(p + p2, N, Q)
	p3_rotated = rotate_homomorphic(p3_rotated, rotation_offset, N) # divide by a redundant scale factor
	p4 = reduce_polynomial((p * p2), N, Q)/scale # divide by a redundant scale factor

	scale_back = scale
	scale = 1
	encoder_1 = BFVEncoder()
	scale = scale_back
	print('unscaled polynomial 1: ' + str(encoder_1.encode(z)))
	print('unscaled polynomial 2: ' + str(encoder_1.encode(z2)))
	#print('unscaled polynomial 1 + 2: ' + str(reduce_polynomial(encoder_1.encode(z) + encoder_1.encode(z2), N, Q)))
	print()
	print('scaled polynomial 1: ' + str(p))
	print('scaled polynomial 2: ' + str(p2))
	print()
	#print("before polynomial reduction of 1 + 2: " + str(p + p2))
	print('after q-reduction of polynomial 1 + polynomial 2: ' + str(reduce_polynomial(p + p2, N, Q)))
	print('after polynomial rotation by ' + str(rotation_offset) + ' positions and q-reduction of polynomial 1 + polynomial 2: ' + str(p3_rotated))
	print()
	#print('before polynomial reduction of 1 * 2: ' + str(p * p2))
	#print('after polynomial reduction of 1 * 2: ' + str(reduce_polynomial(p * p2, N, Q)))
	print('after polynomial reduction & scale normalization of 1 * 2 scale normalization: ' + str(p4))
	
	print()

	decoded_z = encoder.decode(p)
	decoded_z2 = encoder.decode(p2)
	decoded_z3 = encoder.decode(p3)
	decoded_z4 = encoder.decode(p4)
	decoded_z3_rotated = encoder.decode(p3_rotated)
	print("<The decoded BFV plaintext vector>")
	print('decoded vector 1: ' + str(decoded_z))
	print('decoded vector 2: ' + str(decoded_z2))
	print('decoded vector 1 + 2: ' + str(decoded_z3))
	print('decodec vector 1 * 2: ' + str(decoded_z4))
	print('rotated decodec vector 1 + 2 by ' + str(rotation_offset) + ' positions: ' + str(decoded_z3_rotated))
	print()

	if not (z3 == decoded_z3).all(): 
		print("vector 1 + 2 is decoded WRONG!")
		print(z3)
		print(decoded_z3)
		sys.exit(0)
	else:
		print("[FINAL DECODED RESULT " + str(test_index) + "] : " + str(z3) + " == " + str(decoded_z3) + " <---- CORRECT")
		print()
	for i in range(len(z3)//2):
		#print(f'Compare {i} vs {i + rotation_offset % (len(z3)//2)}');
		original_index = i
		rotated_index = (i - rotation_offset) % (len(z3)//2)
		if z3[original_index] != decoded_z3_rotated[rotated_index] or z3[original_index + len(z3)//2] != decoded_z3_rotated[rotated_index + len(z3)//2]:
			print("Rotation Wrong!")
			sys.exit(0)
	print("Rotation is CORRECT")

if __name__ == "__main__":
	test_index = 0
	while True:
		main()
		test_index += 1
		if not is_random_param or test_index == test_count:
			break

	if is_random_param:
		print("Total " + str(test_index) + " Tests Passed")
