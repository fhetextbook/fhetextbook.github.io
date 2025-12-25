import numpy as np
from fastcore.foundation import patch_to
from os import sys
import random
import math
from sympy import nextprime, mod_inverse, primitive_root
import copy
from tfhe_lib import *
import argparse

# We set N and Q as hyper-parameters, based on which P and scale will be automatically configured

Q = 2**30
P = 2**8
K=2**5

#############################

# Test Option 1: Targeted Tests

m1 = 10
m2 = 2

#############################

# Test Option 2: Random Test Params

test_count = 1000

##############################################

np.set_printoptions(suppress=True)


def build_parser():
	 """Build argument parser"""
	 parser = argparse.ArgumentParser(description='Python-based TFHE Demo Library')
	 parser.add_argument('--encrypt', action='store_true', help='Encrytion/decryption test')
	 parser.add_argument('--add-cipher-plain', action='store_true', help='Cipher-plain addition test')
	 parser.add_argument('--add-cipher-cipher', action='store_true', help='Cipher-cipher addition test')
	 parser.add_argument('--mult-cipher-plain', action='store_true', help='Cipher-plain multiplication test')
	 parser.add_argument('--mult-cipher-cipher', action='store_true', help='Cipher-cipher multiplication test')
	 parser.add_argument('--keyswitch', action='store_true', help='Key switch test')
	 #parser.add_argument('--bootstrapping', action='store_true', help='Bootstrapping test')
	 parser.add_argument('--random', action='store_true', help='A bulk of random tests')
	 parser.add_argument('--all', action='store_true', help='All test')
	 return parser


def main(test_index, is_encrypt, run_choice, is_random_param):

	global m1
	global m2
	
	print("=========================================")
	print("'" + run_choice + "'" + " Test")

	if Q <= 0 or (Q & (Q - 1)) != 0:
		print("[Error] Q is not a power of 2: " + str(Q))
		sys.exit()
	if P <= 0 or (P & (P - 1)) != 0:
		print("[Error] P is not a power of 2: " + str(P))
		sys.exit()
	if P >= Q:
		print("[Error] Q is not greater than P: " + str(P) + ", " + str(Q))
		sys.exit()

	scale = Q // P

	# print_cyclotomic_polynomial_info(N=N, P=P, is_complex=False)

	if is_random_param:
		m1 = random.randint(0, math.floor(math.sqrt(P - 1)))
		m2 = random.randint(0, math.floor(math.sqrt(P - 1)))

	if m1 >= P or m2 >= P:
		print('[Error] TFHE\'s message should not overflor P=' + str(P) + ': ' + str(m1) + ", " + str(m2))
		sys.exit()
	
	print("<The TFHE plaintext value>")
	print('message 1: ' + str(m1))
	print('message 2: ' + str(m2))
	print()

	# =========== Encryption Stage ===============

	s = tfhe_key_generate(K)
	a1, b1 = tfhe_encrypt(m1, K, P, Q, s, scale)
	a2, b2 = tfhe_encrypt(m2, K, P, Q, s, scale)
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
		print('m3 = m1 + m2')
		print()
		m3 = (m1 + m2)
		a3, b3 = tfhe_add_cipher_cipher(a1, b1, a2, b2, Q)

	if run_choice == 'add_cipher_plain':
		print("<<Add_Cipher_Plain Test>>")
		print('m3 = m1 + m2')
		print()
		m3 = (m1 + m2)
		a3, b3 = tfhe_add_cipher_plain(a1, b1, m2, Q, scale)

	if run_choice == 'mult_cipher_plain':
		print("<<Multiply_Cipher_Plain Test>>")
		print('m3 = m1 * m2')
		print()
		m3 = (m1 * m2)
		a3, b3 = tfhe_mult_cipher_plain(a1, b1, m2, Q)

	if run_choice == 'keyswitch':
		print("<<Key_Switch Test>>")
		print('m3: m1 encrypted by a new secret key')
		print()
		s_new = tfhe_key_generate(K)
		ksk = tfhe_keyswitch_keygen(s, s_new, K, Q, base=1<<10)
		a3, b3 = tfhe_keyswitch(a1, b1, ksk)
		m3 = m1

	if run_choice == 'mult_cipher_cipher':
		print("<<Multiply_Cipher_Cipher Test>>")
		print('m3: m1 * m2')
		print()
		m3 = m1 * m2

		### <mult-cipher-cipher by the same way as BFV (more efficient)>
		s_kron = np.outer(s, s).flatten() # Flattened tensor product of s
		relin_key = tfhe_keyswitch_keygen(s_kron, s, K, Q)
		a3, b3 = tfhe_mult_cipher_cipher(a1, b1, a2, b2, K, P, Q, scale, relin_key)

		### <mult-cipher-cipher by circuit bootstrapping (less efficient)>
		# bk = tfhe_bootstrap_keygen(s)
		# a3, b3 = tfhe_mult_cipher_cipher_circuit_bootstrapping(a1, b1, a2, b2, K, P, Q, scale, bk)

	if run_choice == 'bootstrapping':
		print("<<Bootstrapping Test>>")
		print('m3: m1 encrypted in a bootstrapped ciphertext')
		print()
		m3 = m1
		a3, b3 = tfhe_bootstrapping(a1, b1, N, P, Q)

	if run_choice != 'encrypt' and (m3 >= P or m3 < 0):
		print('[Error] m3 is out of P: ' + str(m3) + ' vs ' + str(P))
		sys.exit()

	# =========== Decryption Stage ===============
	print('<Decrypted Values>')
	scaling_bits = math.floor(math.log2(scale))
	dec_m1, e1_bits = tfhe_decrypt(a1, b1, P, Q, s, scale, m1)
	dec_m2, e2_bits = tfhe_decrypt(a2, b2, P, Q, s, scale, m2)
	print('- m1:           ', m1)
	print('- decrypted_m1: ', dec_m1)
	print(f'- scaling-factor-bit : noise-bit ratio = {scaling_bits} : {e1_bits}')
	print()
	print('- m2:           ', m2)
	print('- decrypted_m2: ', dec_m2)
	print(f'- scaling-factor-bit : noise-bit ratio = {scaling_bits} : {e2_bits}')
	print()
	if run_choice != 'encrypt':
		if run_choice == 'keyswitch':
			dec_m3, e3_bits = tfhe_decrypt(a3, b3, P, Q, s_new, scale, m3)
		else:
			dec_m3, e3_bits = tfhe_decrypt(a3, b3, P, Q, s, scale, m3)
		print('- m3:           ', m3)
		print('- decrypted_m3: ', dec_m3)
		print(f'- scaling-factor-bit : noise-bit ratio = {scaling_bits} : {e3_bits}')
		print()

	print()
	if m1 != dec_m1 or m2 != dec_m2:
		print("[Error] decoded messages mismatch")
		print("- m1:", m1, dec_m1)
		print("- m2:", m2, dec_m2)
		sys.exit(0)
	if run_choice != 'encrypt' and dec_m3 != m3:
		print("[Error] decoded vectors mismatch")
		print("- m3:", m3, dec_m3)
		sys.exit(0)

	print(f"[TEST {test_index + 1}] '{run_choice}': all values exactly match")
	print()

if __name__ == "__main__":
	args = build_parser().parse_args()
	if not any(vars(args).values()):
		args.random = True
		args.all = True
	if args.all:
		args.random = True
	is_random_param = args.random
	if args.all or args.random and not (args.encrypt or args.add_cipher_cipher or args.add_cipher_plain or args.mult_cipher_cipher or args.mult_cipher_plain or args.keyswitch or args.bootstrapping): 
		args.encrypt = True
		args.add_cipher_cipher = True
		args.add_cipher_plain = True
		args.mult_cipher_cipher = True
		args.mult_cipher_plain = True
		args.keyswitch = True
		args.bootstrapping = True
	run_slot = []
	if args.encrypt:
		run_slot.append('encrypt')
	if args.add_cipher_cipher:
		run_slot.append('add_cipher_cipher')
	if args.add_cipher_plain:
		run_slot.append('add_cipher_plain')
	if args.mult_cipher_cipher:
		run_slot.append('mult_cipher_cipher')
	if args.mult_cipher_plain:
		run_slot.append('mult_cipher_plain')
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
