import numpy as np
from numpy.polynomial import Polynomial
from fastcore.foundation import patch_to
from os import sys
import random
import math
from sympy import nextprime, mod_inverse, primitive_root
import copy
from bgv_lib import *
import argparse
from common import *

# We set N and Q as hyper-parameters, based on which P and scale will be automatically configured

Q_LEVEL = 5
Q_BIT = 30
N = 4
error_tolerance_ratio = 0.01

#############################

# Test Option 1: Targeted Tests

# Testcase 1
bgv_vector1 = np.array([10, 3, 5, 13], dtype=object)
bgv_vector2 = np.array([2, 4, 3, 6], dtype=object)
rotation_offset = 3

# Testcase 2
#N = 8
#bgv_vector1 = np.array([1 + 2j, -3 + 4j, 5 - 6j, -7 -8j], dtype=object)
#bgv_vector2 = np.array([0, 0, 0, 0], dtype=object)
#rotation_offset = 3



#############################

# Test Option 2: Random Test Params

test_count = 1000

##############################################

np.set_printoptions(suppress=True)
is_homogeneous_rotation = True


def build_parser():
	 """Build argument parser"""
	 parser = argparse.ArgumentParser(description='Python-based BGV Demo Library')
	 parser.add_argument('--encode', action='store_true', help='Encoding test')
	 parser.add_argument('--encrypt', action='store_true', help='Encrytion/decryption test')
	 parser.add_argument('--add-cipher-plain', action='store_true', help='Cipher-plain addition test')
	 parser.add_argument('--add-cipher-cipher', action='store_true', help='Cipher-cipher addition test')
	 parser.add_argument('--mult-cipher-plain', action='store_true', help='Cipher-plain multiplication test')
	 parser.add_argument('--mult-cipher-cipher', action='store_true', help='Cipher-cipher multiplication test')
	 parser.add_argument('--keyswitch', action='store_true', help='Key switch test')
	 parser.add_argument('--rotate', action='store_true', help='Rotation test')
	 #parser.add_argument('--bootstrapping', action='store_true', help='Bootstrapping test')
	 parser.add_argument('--random', action='store_true', help='A bulk of random tests')
	 parser.add_argument('--all', action='store_true', help='All test')
	 return parser


def main(test_index, is_encrypt, run_choice, is_random_param):

	global bgv_vector1
	global bgv_vector2

	if run_choice == '':
		run_choice = 'encrypt' if args.encrypt else 'encode'
	
	P = init_P(N=N)
	Q_BASE = (pick_ntt_primes(bits=Q_BIT, num=Q_LEVEL, p_mod=P, max_span=1<<11))
	# print_cyclotomic_polynomial_info(N=N, P=P, is_complex=False)

	budget_bits = math.floor(math.log2(math.prod(Q_BASE) / P))

	scale = P
	Q_BASE_ORIGIN = Q_BASE.copy()
	is_prime.cache = {}
	is_prime.cache[math.prod(Q_BASE)] = False

	if is_random_param:
		rotation_offset = random.randint(-N//2 + 1, N//2 - 1)
		bgv_vector1 = np.empty(N, dtype=object)
		bgv_vector2 = np.empty(N, dtype=object)
		for i in range(N):
			bgv_vector1[i] = centered_mod(random.randint(0, P - 1), P)
			bgv_vector2[i] = centered_mod(random.randint(0, P - 1), P)
	else:
		if len(bgv_vector1) != N or len(bgv_vector2) != N:
			print("[Error] The dimensions of the input vectors are not " + str(N) + ", but " + str(len(bgv_vector1)) + ", " + str(len(bgv_vector2)))
			sys.exit()


	encoder = Encoder(N=N, P=P, scale=scale)

	print("<The BGV plaintext vector to encode>")

	bgv_vector1 = centered_mod_arr(bgv_vector1, P)
	bgv_vector2 = centered_mod_arr(bgv_vector2, P)

	z1 = np.array(bgv_vector1, dtype=object)
	z2 = np.array(bgv_vector2, dtype=object)
	print('vector 1: ' + str(z1))
	print('vector 2: ' + str(z2))
	print()
	# =========== Encryption Stage ===============

	print("<The encoded BGV plaintext polynomial>")
	p1 = encoder.encode(z1)
	p2 = encoder.encode(z2)
	print('encoded polynomial 1: ' + str(p1))
	print('encoded polynomial 2: ' + str(p2))
	print()

	if args.encrypt:
		s = key_generate(N)
		a1, b1 = fhe_encrypt(p1, N, P, Q_BASE, s, encoder.scale)
		a2, b2 = fhe_encrypt(p2, N, P, Q_BASE, s, encoder.scale)
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
		a3, b3 = fhe_add_cipher_cipher(a1, b1, a2, b2, N, Q_BASE)

	if run_choice == 'add_cipher_plain':
		print("<<Add_Cipher_Plain Test>>")
		print('z3 = z1 + z2')
		print()
		z3 = centered_mod_arr(z1 + z2, P)
		p3 = poly_add(p1, p2, N, P)
		a3, b3 = fhe_add_cipher_plain(a1, b1, p2, N, Q_BASE, encoder.scale)

	if run_choice == 'mult_cipher_plain':
		print("<<Multiply_Cipher_Plain Test>>")
		print('z3 = z1 * z2')
		print()
		z3 = centered_mod_arr(z1 * z2, P)
		p3 = poly_mult(p1, p2, N)
		a3, b3 = fhe_mult_cipher_plain(a1, b1, p2, N, Q_BASE)

	if run_choice == 'keyswitch':
		print("<<Key_Switch Test>>")
		print('z3: z1 encrypted by a new secret key')
		print()
		s_new = key_generate(N)
		ksk = fhe_keyswitch_keygen(s, s_new, N, math.prod(Q_BASE), base=1<<10, bgv_scale=P)
		a3, b3 = fhe_keyswitch(a1, b1, ksk)
		p3 = p1
		z3 = z1.copy()

	if run_choice == 'mult_cipher_cipher':
		print("<<Multiply_Cipher_Cipher Test>>")
		print('z3: z1 * z2')
		print()
		z3 = centered_mod_arr(z1 * z2, P)
		p3 = poly_mult(p1, p2, N)

		s2 = poly_mult(s, s, N, math.prod(Q_BASE))
		relin_key = fhe_keyswitch_keygen(s2, s, N, math.prod(Q_BASE), bgv_scale=P)
		a3, b3 = fhe_mult_cipher_cipher(a1, b1, a2, b2, N, P, Q_BASE, encoder.scale, relin_key)

	if run_choice == 'rotate':
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
		p3 = poly_rotate(p1, math.prod(Q_BASE_ORIGIN), N, rotation_offset) # divide by a redundant scale factor

		rotate_ksk = poly_mod(poly_rotate(s, math.prod(Q_BASE), N, rotation_offset), N, math.prod(Q_BASE))
		galois_rotation_key = fhe_keyswitch_keygen(rotate_ksk, s, N, math.prod(Q_BASE), base=1 << 10, bgv_scale=P)
		a3, b3 = fhe_rotate(a1, b1, N, Q_BASE, rotation_offset, galois_rotation_key)

	if run_choice == 'bootstrapping':
		print("<<Bootstrapping Test>>")
		print('z3: z1 encrypted in a bootstrapped ciphertext')
		print()
		z3 = z1.copy()
		p3 = p1.copy()
		a3, b3 = fhe_bootstrapping(a1, b1, N, Q_BASE)

	# =========== Decryption Stage ===============
	if args.encrypt:
		print('<Decrypted Polynomials>')
		#scaling_bits = math.floor(math.log2(scale))
		dec_p1, e1_bits = fhe_decrypt(a1, b1, N, P, Q_BASE_ORIGIN, s, encoder.scale, p1)
		dec_p2, e2_bits = fhe_decrypt(a2, b2, N, P, Q_BASE_ORIGIN, s, encoder.scale, p2)
		print('- p1:			  ', p1)
		print('- decrypted_p1: ', dec_p1)
		print(f'- Q-budget-bits : maximum noise-bits ratio = {budget_bits} : {e1_bits}')
		print()
		print('- p2:			  ', p2)
		print('- decrypted_p2: ', dec_p2)
		print(f'- Q-budget-bits : maximum noise-bits ratio = {budget_bits} : {e2_bits}')
		print()
		if run_choice != 'encode' and run_choice != 'encrypt':
			if run_choice == 'keyswitch':
				dec_p3, e3_bits = fhe_decrypt(a3, b3, N, P, Q_BASE, s_new, encoder.scale, p3)
			else:
				dec_p3, e3_bits = fhe_decrypt(a3, b3, N, P, Q_BASE, s, encoder.scale, p3)
			print('- p3:			  ', p3)
			print('- decrypted_p3: ', dec_p3)
			print(f'- Q-budget-bits : noise-bits ratio = {budget_bits} : {e3_bits}')
			print()
	else:
		dec_p1 = p1
		dec_p2 = p2

	decoded_z1 = encoder.decode(dec_p1)
	decoded_z2 = encoder.decode(dec_p2)
	
	print('<Decoded vectors>')
	print('- z1:			', z1)
	print('- decoded_z1: ', decoded_z1)
	print()
	print('- z2:			', z2)
	print('- decoded_z2: ', decoded_z2)
	print()
	if run_choice != 'encode' and run_choice != 'encrypt':
		decoded_z3 = encoder.decode(dec_p3)
		print('- z3:			', z3)
		print('- decoded_z3: ', decoded_z3)
		print()
	print()
	if (not np.array_equal(z1, decoded_z1)) or (not np.array_equal(z2, decoded_z2)):
		print("[Error] decoded vectors mismatch")
		print("- z1:", z1, decoded_z1)
		print("- z2:", z2, decoded_z2)
		sys.exit(0)
	if run_choice != 'encode' and run_choice != 'encrypt' and not np.array_equal(z3, decoded_z3):
		print("[Error] decoded vectors mismatch")
		print("- z3:", z3, decoded_z3)
		sys.exit(0)

	print(f"[TEST {test_index + 1}] '{run_choice}': all values exactly match")
	print()
	return



if __name__ == "__main__":
	args = build_parser().parse_args()
	if not any(vars(args).values()):
		args.random = True
		args.all = True
	if args.all:
		args.random = True
	is_random_param = args.random
	if args.all or args.random and not (args.encode or args.encrypt or args.add_cipher_cipher or args.add_cipher_plain or args.mult_cipher_cipher or args.mult_cipher_plain or args.keyswitch or args.rotate): 
		args.encrypt = True
		args.add_cipher_cipher = True
		args.add_cipher_plain = True
		args.mult_cipher_cipher = True
		args.mult_cipher_plain = True
		args.rotate = True
		args.keyswitch = True
	if args.encrypt:
		args.encode = True
	if args.add_cipher_cipher or args.add_cipher_plain or args.mult_cipher_cipher or args.mult_cipher_plain or args.keyswitch or args.rotate:
		args.encrypt = True
	run_slot = []
	if args.add_cipher_cipher:
		run_slot.append('add_cipher_cipher')
	if args.add_cipher_plain:
		run_slot.append('add_cipher_plain')
	if args.mult_cipher_cipher:
		run_slot.append('mult_cipher_cipher')
	if args.mult_cipher_plain:
		run_slot.append('mult_cipher_plain')
	if args.rotate:
		run_slot.append('rotate')
	if args.keyswitch:
		run_slot.append('keyswitch')
	#if args.bootstrapping:
	#	run_slot = ['bootstrapping']

	test_index = 0
	while True:
		main(test_index, args.encrypt, random.choice(run_slot) if run_slot else '', is_random_param)
		test_index += 1
		if not is_random_param or test_index == test_count:
			break

	if is_random_param:
		print("Total " + str(test_index) + " Tests Passed")
