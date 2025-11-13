from sympy import randprime
import random
import sys
import numpy as np
import math
from sympy import nextprime, mod_inverse, primitive_root
from functools import reduce
import operator
from collections import Counter
np.set_printoptions(threshold=np.inf)

CIPHER_BITS = 14
SM_MODULUS_BITS = 4
EX_MODULUS_BITS = 10
PLAIN_BITS = 8
MODULUS_COUNT = 5
N = 1 << 11
prime_e = 3

def centered_mod(x, q):
    """
    Centered residue of x (scalar or array‑like) modulo q.

    Returned range:
        (-⌊q/2⌋ , ⌈q/2⌉]   for any positive integer q.
    """

    if q <= 0:
        raise ValueError("Modulus q must be positive.")

    # ---------- scalar fast‑path ----------
    if np.isscalar(x):
        r = x % q                     # bring into [0, q-1]
        if r > q // 2:                # fold the upper half down
            r -= q
        return int(r)

    # ---------- vectorised path ----------
    arr = np.asarray(x)               # 1‑D, 2‑D, … anything works
    r = np.mod(arr, q)                # same shape as arr
    mask = r > q // 2
    r = r.astype(object, copy=False)  # allow big ints if q is big
    r[mask] -= q
    return r if isinstance(x, np.ndarray) else type(x)(r.tolist())

def random_prime_sympy(bits):
	 lower_bound = 1 << (bits - 1)
	 upper_bound = (1 << bits) - 1
	 return randprime(lower_bound, upper_bound)

def find_duplicates(arr):
	 return [item for item, count in Counter(arr).items() if count > 1]



def mod_inverse(x, q):
    if q <= 0:
        raise ValueError("Modulus q must be positive.")

    # ---- scalar fast path -------------------------------------------------
    if np.isscalar(x):
        return _scalar_inv_centered(int(x), q)

    # ---- vectorised path --------------------------------------------------
    arr  = np.asarray(x, dtype=object)          # keeps big ints intact
    invs = np.empty_like(arr)

    it = np.nditer(arr, flags=["multi_index"])
    while not it.finished:
        invs[it.multi_index] = _scalar_inv_centered(int(it[0]), q)
        it.iternext()

    return invs if isinstance(x, np.ndarray) else type(x)(invs.tolist())


def _scalar_inv_centered(a, q):
    a_mod = a % q
    if math.gcd(a_mod, q) != 1:
        raise ValueError(f"{a} has no inverse modulo {q}")

    # Python 3.8+:  pow(base, -1, mod) gives modular inverse
    inv = pow(a_mod, -1, q)
    return centered_mod(inv, q)



def poly_reduce(prod, modulus):
	# Wrap around for x^N - 1
	for i in range(N, len(prod)):
		prod[i - N] += prod[i]
	reduced = prod[:N]		  # shape (16384,)
	for i in range(0, N):
		prod[i] = centered_mod(prod[i], modulus)
	return reduced

def poly_round(poly):
	for i in range(0, N):
		poly[i] = round(poly[i])
	return poly

def rns_recover(M_rns, q_base):
	M_recovered = 0
	q = reduce(operator.mul, q_base, 1)
	multiple=0	
	for i in range(0, len(q_base)):
		y = (q // q_base[i])
		z = mod_inverse(centered_mod(y, q_base[i]), q_base[i])
		multiple += y * z
		M_recovered += centered_mod((centered_mod(M_rns[i] * z, q_base[i])) * y, q)
	return centered_mod(M_recovered, q)

def FastBConv(val_rns, q_base, b_base, is_exact=False):
	merge = 0
	q = reduce(operator.mul, q_base, 1)
	if len(val_rns) != len(q_base):
		print(f'Error: val_rns={len(val_rns)} is not the same length as q_base={len(q_base)}')
		raise Exception
	# for each modulus's polynomial
	for i in range(0, len(val_rns)):
		y = int(q // q_base[i])
		z = mod_inverse(y, q_base[i])  
		merge += (centered_mod(val_rns[i] * z, q_base[i])) * y 

	if is_exact:
		merge = centered_mod(merge, q)

	converted_val_rns = []	
	for i in range(0, len(b_base)):
		converted_val_rns.append(centered_mod(merge, b_base[i]))
	'''
	merge_test = centered_mod(merge, q)
	print("DIFF: ", (merge - merge_test) / q)
	merge2 = rns_recover(converted_val_rns, b_base)
	b = reduce(operator.mul, b_base, 1)
	print("DIFF 2: ", (merge2 - (centered_mod(rns_recover(val_rns, q_base), b))) / q)
	sys.exit()	
	print("Check Diff :", centered_mod((rns_recover(converted_poly_rns, b_base[len(q_base):]) - (centered_mod(rns_recover(poly_rns, q_base), b))), q))
	'''
	return converted_val_rns;

def ModRaiseRNS(poly_rns, q_base, b_base):
	converted_poly_rns = FastBConv(poly_rns, q_base, b_base[len(q_base):])
	merged_poly_rns = poly_rns + converted_poly_rns
	'''
	print(q_base)
	print(b_base[len(q_base):])
	b = reduce(operator.mul, b_base, 1)
	print("Check Diff :", (rns_recover(converted_poly_rns, b_base[len(q_base):]) - (centered_mod(rns_recover(poly_rns, q_base), b))) / q)
	print("Check Diff 2 :", (rns_recover(merged_poly_rns, b_base) - (centered_mod(rns_recover(poly_rns, q_base), b))) / q)
	sys.exit()
	'''
	return merged_poly_rns 


def ModDropRNS(poly_rns, q_base, subset_base):
	converted_poly_rns = FastBConv(poly_rns, q_base, subset_base)
	return converted_poly_rns

def ModSwitchRNS(poly_rns, first_second_base, first_base):
	first_poly_rns = poly_rns[:len(first_base)]
	second_poly_rns = poly_rns[len(first_base):]
	second_base = first_second_base[len(first_base):]
	converted_second_poly_rns = FastBConv(second_poly_rns, second_base, first_base)
	switched_poly_rns = []
	b = reduce(operator.mul, b_base, 1)
	for i in range(0, len(first_base)):
		switched_poly_rns.append(centered_mod((mod_inverse(b, first_base[i]) * (first_poly_rns[i] - converted_second_poly_rns[i])), first_base[i]))
	return switched_poly_rns

def SmallMont(poly_rns, augmented_b_base, sm_modulus, q_base):
	recovered_poly_rns = rns_recover(poly_rns, augmented_b_base)
	q = reduce(operator.mul, q_base, 1)
	augmented_b = reduce(operator.mul, augmented_b_base, 1)
	b = augmented_b // sm_modulus
	final_sm_reduced_poly_rns = []
	common_poly = centered_mod(recovered_poly_rns * mod_inverse(q, sm_modulus), sm_modulus)
	sm_reduced_poly_rns_debug = centered_mod((recovered_poly_rns - (common_poly * q)) // sm_modulus, b)
	b_base = []
	for i in range(0, len(augmented_b_base)):
		if augmented_b_base[i] == sm_modulus:
			continue
		b_base.append(augmented_b_base[i])
		reduced_recovered_poly_rns = centered_mod(recovered_poly_rns, augmented_b_base[i])
		for j in range(0, N):
			if reduced_recovered_poly_rns[j] != poly_rns[i][j]:
				print("Error: reduced poly_rns is not the same as manually computed one:")
				sys.exit()
		reduced_q = centered_mod(q, augmented_b_base[i])
		final_sm_reduced_poly_rns.append(centered_mod(((centered_mod(recovered_poly_rns, augmented_b_base[i]) - reduced_q * common_poly) * mod_inverse(sm_modulus, augmented_b_base[i])), augmented_b_base[i]))

	recovered_final_sm_reduced_poly = rns_recover(final_sm_reduced_poly_rns, b_base)

	for i in range(0, N):
		if recovered_final_sm_reduced_poly[i] != sm_reduced_poly_rns_debug[i]:
			print(f'Error: sm_reduced poly does not match')
			print(recovered_final_sm_reduced_poly[i])
			print(sm_reduced_poly_rns_debug[i])
			sys.exit()
	return final_sm_reduced_poly_rns

def FastBConvEx(poly_rns, b_base, ex_modulus, q_base):
	x_hat_poly_rns = ModDropRNS(poly_rns, b_base + [ex_modulus], b_base)
	x_alpha_poly_rns = ModDropRNS(poly_rns, b_base + [ex_modulus], ex_modulus)
	b = reduce(operator.mult, b_base)
	gamma = centered_mod((FastBConv(x_hat_poly_rns, b_base, ex_modulus) - x_alpha_poly_rns) * mod_inverse(b, ex_modulus), ex_modulus)
	converted_x_hat_poly_rns = FastBConv(x_hat_poly_rns, b_base, q_base)
	
	converted_exact_poly_rns = []
	for i in range(0, len(q_base)):
		b_converted = centered_mod(b, q_base[i])
		converted_exact_poly_rns.append(centered_mod((converted_x_hat_poly_rns[i] - gamma * b_converted) , q_base[i]))
	return converted_exact_poly_rns

all_moduli = []

# Find unique ciphertext moduli
while True:
	q_base=[int(random_prime_sympy(CIPHER_BITS)) for _ in range(0, MODULUS_COUNT)]
	if len(find_duplicates(q_base)) == 0:
		all_moduli += q_base
		break

# Find unique axilliary FastBConv moduli
while True:
	b_base=[int(random_prime_sympy(CIPHER_BITS + 1)) for _ in range(0, MODULUS_COUNT)]
	b_sk = int(random_prime_sympy(CIPHER_BITS + 1))
	if len(find_duplicates(b_base + [b_sk] + all_moduli)) == 0:
		all_moduli += b_base + [b_sk]
		break

# Find a unique axilliary SmallMont modulus
while True:
	sm_modulus = int(random_prime_sympy(SM_MODULUS_BITS))
	if len(find_duplicates([sm_modulus] + all_moduli)) == 0:
		all_moduli += [sm_modulus]
		break

# Find a unique axilliary FastBConvEx modulus
while True:
	ex_modulus = int(random_prime_sympy(EX_MODULUS_BITS))
	if len(find_duplicates([ex_modulus] + all_moduli)) == 0:
		all_moduli += [ex_modulus]
		break

q = reduce(operator.mul, q_base, 1)
b = reduce(operator.mul, b_base, 1)
p = random_prime_sympy(PLAIN_BITS)
scale = math.floor(q / p)


print('- ciphertext modulus base: ' + str(q_base))
print('- plaintext modulus: ' + str(p))
print('- sm modulus: ' + str(sm_modulus))
print('- ex modulus: ' + str(ex_modulus))

A = np.array([random.randrange(q) for _ in range(N)], dtype=object)
M = np.array([random.randrange(p) for _ in range(N)], dtype=object)
E = np.array([random.randrange(int(p / 4)) for _ in range(N)], dtype=object)
S = np.array([random.randrange(3) for _ in range(N)], dtype=object)

S -= 1

A_rns = []
E_rns = []
M_rns = []
S_rns = []
B_rns = []

print("M[0]: " + str(M[0]))
for i in range(0, MODULUS_COUNT):
	A_rns.append(centered_mod(A, q_base[i]))
	E_rns.append(centered_mod(E, q_base[i]))
	M_rns.append(centered_mod(M, q_base[i]))
	S_rns.append(centered_mod(S, q_base[i]))
	print(f'M_rns[{i}][0]: {M_rns[i][0]}')

### Encrypt Big
B = poly_reduce(scale * M + E - poly_reduce(np.convolve(A, S), q), q)

### Encrypt RNS
for i in range(0, MODULUS_COUNT):
	temp = poly_reduce(np.convolve(A_rns[i], S_rns[i]), q_base[i])
	B_rns.append(poly_reduce(scale * M_rns[i] + E_rns[i] - temp, q_base[i])) 

### Decrypt Big
decrypted_M = centered_mod(poly_round(poly_reduce((poly_reduce(np.convolve(A, S), q) + B), q) / scale), p)
print("Decrypted M[0]: " + str(decrypted_M[0]))

### Decrypt RNS & Recover
decrypted_M_rns = []
for i in range(0, MODULUS_COUNT):
	temp = poly_reduce(np.convolve(A_rns[i], S_rns[i]), q_base[i])
	decrypted_M_rns.append(temp + B_rns[i])
decrypted_M_rns_recovered = centered_mod(poly_round(rns_recover(decrypted_M_rns, q_base) / scale), p)
print("Decrypted M_rns_recovered[0]: " + str(decrypted_M_rns_recovered[0]))


recovered_A = rns_recover(A_rns, q_base)

converted_A_rns = FastBConv(A_rns, q_base, b_base)
converted_A_rns_debug = FastBConv(A_rns, q_base, b_base, True)

converted_A_rns_recovered = rns_recover(converted_A_rns, b_base)
converted_A_rns_debug_recovered = rns_recover(converted_A_rns_debug, b_base)

raised_A_rns = ModRaiseRNS(A_rns, q_base, q_base + b_base)
recovered_raised_A = rns_recover(raised_A_rns, q_base + b_base)

raised_switched_A_rns = ModSwitchRNS(raised_A_rns, q_base + b_base, q_base)
recovered_raised_switched_A = rns_recover(raised_switched_A_rns, q_base)

multiplied_A_rns = []
for i in range(0, len(q_base)):
	multiplied_A_rns.append(centered_mod(A_rns[i] * sm_modulus, q_base[i]))

multiplied_converted_A_rns = FastBConv(multiplied_A_rns, q_base, b_base + [sm_modulus])
A_sm_reduced_rns = SmallMont(multiplied_converted_A_rns, b_base + [sm_modulus], sm_modulus, q_base)

# Verify that the SmallMontgomery reduction's q-overflow is only {-1, 0, 1}
sm_q_overflow_list = centered_mod(centered_mod(rns_recover(A_rns, q_base) - rns_recover(A_sm_reduced_rns, b_base), b) * mod_inverse(q, b), b)
for i in sm_q_overflow_list:
	if i > 1 or i < -1:
		print("Error: SmallMontgomery reduction's q-overflow is larger than abs(1): " + str(i))
		sys.exit()
print(sm_q_overflow_list)
print("Small Montgomery reduction's all q-overflows are either {-1, 0, 1}")


