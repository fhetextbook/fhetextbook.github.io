import numpy as np
from numpy.polynomial import Polynomial
from fastcore.foundation import patch_to
from os import sys
import random
import math
from sympy import nextprime, mod_inverse, primitive_root
import copy
from bfv_lib import *
import argparse

# First we set the parameters


#############################

# Test Option 1: Targeted Test Params

N = 4
P = 17
Q = 2**30

bfv_vector1 = np.array([10, 3, 5, 13], dtype=object)
bfv_vector2 = np.array([2, 4, 3, 6], dtype=object)

rotation_offset = 3

#############################

# Test Option 2: Random Test Params

test_count = 1000

N_LOWER = 8
N_UPPER = 8

P_LOWER = 35
P_UPPER = 10000

test_index = 0

##############################################

np.set_printoptions(suppress=True)

is_homogeneous_rotation = True

def InitParam(bfv_vector1, bfv_vector2, is_random_param):
	global N
	global P
	global Q
	global scale
	M = N*2
	global N_inverse
	global primitive_root_of_unity
	#global scale_inverse

	if is_random_param:
		# Decide the exponent N
		while True:
			power = random.randint(N_LOWER, N_UPPER)
			if IsPowerOf2(power):
				 N = power
				 P = random.randint(P_LOWER, P_UPPER)
				 M = N * 2
				 break
	else:
		if len(bfv_vector1) != N or len(bfv_vector2) != N:
		  print("Error: The dimensions of the input vectors are not " + str(N) + ", but " + str(len(bfv_vector1)) + ", " + str(len(bfv_vector2)))
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

	N_inverse = ModInv(N, P) # If P is a prime, the multiplicative inverse of N (< P) is guaranteed to exist (by Fermit's Little theorem)
	scale = math.floor(Q / P)
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
					root_values.append(i)	# Then, i is the primitive M-th root of unity
		#if not is_random_param: 
		if len(root_values) == 0 and not is_random_param: 
			print("Modulus P=" + str(P) + " has no primitive " + str(M) + "-th root of unity")
			sys.exit(0)

	primitive_root_of_unity = min(root_values)
	'''

#  '''
	### [OPTION 2]: An alternative efficient way of computing the primitive (M=2N)-th root of unity
	mult_group_generator = primitive_root(P)  # prim_root is a number such that Order(prim_root) = totient(p) mod p
	primitive_root_of_unity = pow(mult_group_generator, (P - 1) // M, P) # The (M=2N)-th primitive root of unity (i.e., solution for both X^N + 1 = 0 mod p AND X^{2N} = 1 mod P)
	next_root = primitive_root_of_unity
	while True:
		root_values.append(next_root)
		next_root = next_root * (primitive_root_of_unity ** 2) % P
		if next_root == primitive_root_of_unity:
			break
#  '''

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
		bfv_vector1 = np.empty(N, dtype=object)
		bfv_vector2 = np.empty(N, dtype=object)
		for i in range(N):
			bfv_vector1[i] = centered_mod(random.randint(0, P - 1), P)
			bfv_vector2[i] = centered_mod(random.randint(0, P - 1), P)
	return (bfv_vector1,  bfv_vector2)

def build_parser():
	 """Build argument parser"""
	 parser = argparse.ArgumentParser(description='Python-based BFV Demo Library')
	 parser.add_argument('--encoding', action='store_true', help='Encoding test')
	 parser.add_argument('--encrypt', action='store_true', help='Encrytion/decryption test')
	 parser.add_argument('--add-cipher-plain', action='store_true', help='Cipher-plain addition test')
	 parser.add_argument('--add-cipher-cipher', action='store_true', help='Cipher-cipher addition test')
	 parser.add_argument('--mult-cipher-plain', action='store_true', help='Cipher-plain multiplication test')
	 parser.add_argument('--mult-cipher-cipher', action='store_true', help='Cipher-cipher multiplication test')
	 parser.add_argument('--keyswitch', action='store_true', help='Key switch test')
	 parser.add_argument('--rotate', action='store_true', help='Rotation test')
	 parser.add_argument('--bootstrapping', action='store_true', help='Bootstrapping test')
	 parser.add_argument('--random', action='store_true', help='A bulk of random tests')
	 parser.add_argument('--all', action='store_true', help='All test')
	 return parser


def main(is_encrypt, run_choice, is_random_param):

	global M
	global primitive_root_of_unity
	global bfv_vector1
	global bfv_vector2
	bfv_vector1, bfv_vector2 = InitParam(bfv_vector1, bfv_vector2, is_random_param)

	encoder = BFVEncoder(P, Q, N, is_homogeneous_rotation)

	bfv_vector1 = centered_mod_arr(bfv_vector1, P)
	bfv_vector2 = centered_mod_arr(bfv_vector2, P)

	print("<The BFV plaintext vector to encode>")
	z1 = np.array(bfv_vector1, dtype=object)
	z2 = np.array(bfv_vector2, dtype=object)
	print('vector 1: ' + str(z1))
	print('vector 2: ' + str(z2))
	print()
	# =========== Encryption Stage ===============

	print("<The encoded BFV plaintext polynomial>")
	p1 = encoder.encode(z1)
	p2 = encoder.encode(z2)
	print('encoded polynomial 1: ' + str(p1))
	print('encoded polynomial 2: ' + str(p2))
	print()

	if args.encrypt:
		s = key_generate(N)
		a1, b1 = fhe_encrypt(p1, N, P, Q, s, encoder.scale)
		dec_p1, e_bits = fhe_decrypt(a1, b1, N, P, Q, s, scale, p1)
		a2, b2 = fhe_encrypt(p2, N, P, Q, s, encoder.scale)
		dec_p2, e_bits = fhe_decrypt(a2, b2, N, P, Q, s, encoder.scale, p2)
		print('<Encrypted ciphertexts>')
		print('A1: ', a1)
		print('B1: ', a2)
		print()
		print('A2: ', b1)
		print('B2: ', b2)
		print()

	# =========== Computation Stage ===============

	if run_choice == 'add_cipher_cipher':
		print("<<Add_Cipher_Cipher Test>>")
		print('z3 = z1 + z2')
		print()
		z3 = centered_mod_arr(z1 + z2, P)
		p3 = poly_add(p1, p2, N, P)
		a3, b3 = fhe_add_cipher_cipher(a1, b1, a2, b2, N, Q)

	if run_choice == 'add_cipher_plain':
		print("<<Add_Cipher_Plain Test>>")
		print('z3 = z1 + z2')
		print()
		z3 = centered_mod_arr(z1 + z2, P)
		p3 = poly_add(p1, p2, N, P)
		a3, b3 = fhe_add_cipher_plain(a1, b1, p2, N, Q, encoder.scale)

	if run_choice == 'mult_cipher_plain':
		print("<<Multiply_Cipher_Plain Test>>")
		print('z3 = z1 * z2')
		print()
		z3 = centered_mod_arr(z1 * z2, P)
		p3 = poly_mult(p1, p2, N, P)
		a3, b3 = fhe_mult_cipher_plain(a1, b1, p2, N, Q)


	if run_choice == "rotate":
		print("<<Rotation Test>>")
		rotation_offset = random.randint(0, N // 2 - 1)
		print('z3: rotate z1 to the left by ' + str(rotation_offset) + ' slots left-half and right-half, each')
		print()
		z3 = z1.copy()
		for i in range(len(z1)//2):
			original_index = i
			rotated_index = (i - rotation_offset) % (len(z1)//2)
			z3[rotated_index] = z1[original_index] 
			z3[rotated_index + len(z3)//2] = z1[original_index + len(z3)//2]
		p3 = poly_rotate(p1, Q, N, rotation_offset) # divide by a redundant scale factor

		rotate_ksk = poly_mod(poly_rotate(s, Q, N, rotation_offset), N, Q)
		galois_rotation_key = fhe_keyswitch_keygen(rotate_ksk, s, N, Q, base=1 << 10)
		a3, b3 = fhe_rotate(a1, b1, N, Q, rotation_offset, galois_rotation_key)

	if run_choice == 'keyswitch':
		print("<<Key_Switch Test>>")
		print('z3: z1 encrypted by a new secret key')
		print()
		s_new = key_generate(N)
		ksk = fhe_keyswitch_keygen(s, s_new, N, Q, base=1<<10)
		a3, b3 = fhe_keyswitch(a1, b1, ksk)
		p3 = p1
		z3 = z1.copy()

	if run_choice == 'mult_cipher_cipher':
		print("<<Multiply_Cipher_Cipher Test>>")
		print('z3: z1 * z2')
		print()
		z3 = centered_mod_arr(z1 * z2, P)
		p3 = poly_mult(p1, p2, N, P)

		relin_key = fhe_relin_keygen(s, N, Q)
		a3, b3 = fhe_mult_cipher_cipher(a1, b1, a2, b2, N, P, Q, encoder.scale, relin_key)

	if run_choice == 'bootstrapping':
		print("<<Bootstrapping Test>>")
		print('z3: z1 encrypted in a bootstrapped ciphertext')
		print()
		z3 = z1.copy()
		p3 = p1.copy()
		a3, b3 = fhe_bootstrapping(a1, b1, N, P, Q)

	# =========== Decryption Stage ===============
	if args.encrypt:
		print('<Decrypted Polynomials>')
		scaling_bits = math.floor(math.log2(scale))
		dec_p1, e1_bits = fhe_decrypt(a1, b1, N, P, Q, s, encoder.scale, p1)
		dec_p2, e2_bits = fhe_decrypt(a2, b2, N, P, Q, s, encoder.scale, p2)
		print('- p1:           ', p1)
		print('- decrypted_p1: ', dec_p1)
		print(f'- scaling_factor_bit / noise_bit ratio = {scaling_bits} : {e1_bits}')
		print()
		print('- p2:           ', p2)
		print('- decrypted_p2: ', dec_p2)
		print(f'- scaling_factor_bit / noise_bit ratio = {scaling_bits} : {e2_bits}')
		print()
		if run_choice != '':
			if run_choice == 'keyswitch':
				dec_p3, e3_bits = fhe_decrypt(a3, b3, N, P, Q, s_new, encoder.scale, p3)
			else:
				dec_p3, e3_bits = fhe_decrypt(a3, b3, N, P, Q, s, encoder.scale, p3)
			print('- p3:           ', p3)
			print('- decrypted_p3: ', dec_p3)
			print(f'- scaling_factor_bit / noise_bit ratio = {scaling_bits} : {e3_bits}')
			print()
	else:
		dec_p1 = p1
		dec_p2 = p2

	decoded_z1 = encoder.decode(dec_p1)
	decoded_z2 = encoder.decode(dec_p2)
	
	print('<Decoded vectors>')
	print('- z1:         ', z1)
	print('- decoded_z1: ', decoded_z1)
	print()
	print('- z2:         ', z2)
	print('- decoded_z2: ', decoded_z2)
	print()
	if run_choice != '':
		decoded_z3 = encoder.decode(dec_p3)
		print('- z3:         ', z3)
		print('- decoded_z3: ', decoded_z3)
		print()

	print()
	if (not np.array_equal(z1, decoded_z1)) or (not np.array_equal(z2, decoded_z2)):
		print("Error: decoded vectors mismatch")
		print("- z1:", z1, decoded_z1)
		print("- z2:", z2, decoded_z2)
		sys.exit(0)
	if run_choice != '' and not np.array_equal(z3, decoded_z3):
		print("Error: decoded vectors mismatch")
		print("- z3:", z3, decoded_z3)
		sys.exit(0)

	return
	if run_choice != '':
		decoded_z3 = encoder.decode(p3)
		if (not np.array_equal(z3, decoded_z3)): 
			print("Error: decoded vectors mismatch")
			print('z3:', z3, decoded_z3)
			sys.exit(0)

	return

	#decoded_z3_rotated = encoder.decode(p3_rotated)

	#print("before polynomial reduction of 1 + 2: " + str(p + p2))
	print('after q-reduction of polynomial 1 + polynomial 2: ' + str(p3))
	print('after polynomial rotation by ' + str(rotation_offset) + ' positions and q-reduction of polynomial 1 + polynomial 2: ' + str(p3_rotated))
	print()	
	print()
	print('decoded vector 1 + 2: ' + str(decoded_z3))
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
	args = build_parser().parse_args()
	if not any(vars(args).values()):
		args.random = True
		args.all = True
	if args.all:
		args.random = True
	is_random_param = args.random
	if args.all or args.random and not (args.encoding or args.encrypt or args.add_cipher_cipher or args.add_cipher_plain or args.mult_cipher_cipher or args.mult_cipher_plain or args.keyswitch or args.bootstrapping or args.rotate): 
		args.encrypt = True
		args.add_cipher_cipher = True
		args.add_cipher_plain = True
		args.mult_cipher_cipher = True
		args.mult_cipher_plain = True
		args.rotate = True
		args.keyswitch = True
		args.bootstrapping = True
	if args.encrypt:
		args.encoding = True
	if args.add_cipher_cipher or args.add_cipher_plain or args.mult_cipher_cipher or args.mult_cipher_plain or args.keyswitch or args.bootstrapping or args.rotate:
		args.encrypt = True
	run_slot = []
	if args.add_cipher_cipher:
		run_slot = ['add_cipher_cipher']
	if args.add_cipher_plain:
		run_slot = ['add_cipher_plain']
	if args.mult_cipher_cipher:
		run_slot = ['mult_cipher_cipher']
	if args.mult_cipher_plain:
		run_slot = ['mult_cipher_plain']
	if args.rotate:
		run_slot = ['rotate']
	if args.keyswitch:
		run_slot = ['keyswitch']
	#if args.bootstrapping:
	#	run_slot = ['bootstrapping']

	test_index = 0
	while True:
		main(args.encrypt, random.choice(run_slot) if run_slot else '', is_random_param)
		test_index += 1
		if not is_random_param or test_index == test_count:
			break

	if is_random_param:
		print("Total " + str(test_index) + " Tests Passed")
