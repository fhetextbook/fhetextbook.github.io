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

N = 1 << 11
PLAIN_BITS = 4
CIPHER_BITS = 11
MODULUS_COUNT = 2
SWITCH_CIPHER_BITS = 10
SWITCH_MODULUS_COUNT = 2

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




def poly_degree_reduce(prod):
	# Wrap around for x^N - 1
	for i in range(N, len(prod)):
		prod[i - N] += prod[i]
	reduced = prod[:N]		  # shape (16384,)
	return reduced

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


def ModSwitchBGV(poly, plain_modulus, cur_modulus, new_modulus):
	poly_divide_round = np.zeros_like(poly)
	for i in range(0, len(poly)):
		poly_divide_round[i] = int(np.round(float(poly[i]) * (float(new_modulus) / float(cur_modulus))))
	epsilon = new_modulus * poly - cur_modulus * poly_divide_round
	H = poly_reduce(epsilon * mod_inverse(cur_modulus, plain_modulus), plain_modulus)
	poly_switch = poly_reduce(poly_divide_round + H, new_modulus)
	return poly_switch

def ModSwitchBGV_debug(poly, plain_modulus, cur_modulus, new_modulus):
	poly_divide_round = np.zeros_like(poly)
	for i in range(0, len(poly)):
		poly_divide_round[i] = int(np.round(float(poly[i]) * (float(new_modulus) / float(cur_modulus))))
	epsilon = new_modulus * poly - cur_modulus * poly_divide_round
	H = poly_reduce(epsilon * mod_inverse(cur_modulus, plain_modulus), plain_modulus)
	poly_switch = poly_reduce(poly_divide_round + H, new_modulus)
	#print(poly_switch)
	#if (centered_mod(poly, plain_modulus) - centered_mod(poly_switch, plain_modulus) != 0).any():
	#		print ("Plain Modulus False")
	#		sys.exit()
		#print(f'INverse: {centered_mod(mod_inverse(prime_factor, new_modulus), plain_modulus)}')
	return (poly_switch, poly_divide_round, H)

def ModSwitchBGV2(poly, plain_modulus, cur_modulus, new_modulus):
	poly_divide_round = poly
	for i in range(0, len(poly)):
		poly_divide_round[i] = int(np.round(float(poly[i]) * (float(new_modulus) / float(cur_modulus))))
	epsilon = new_modulus * poly - cur_modulus * poly_divide_round
	H = poly_reduce(epsilon * mod_inverse(cur_modulus, plain_modulus), plain_modulus)
	poly_switch = poly_divide_round + H
	#print(poly_switch)
	#if (centered_mod(poly, plain_modulus) - centered_mod(poly_switch, plain_modulus) != 0).any():
	#		print ("Plain Modulus False")
	#		sys.exit()
		#print(f'INverse: {centered_mod(mod_inverse(prime_factor, new_modulus), plain_modulus)}')
	return poly_switch


all_moduli = []
q_base = []
q_switch_base = []
p = random_prime_sympy(PLAIN_BITS)

# Find unique ciphertext moduli such that each moduli is 1 mod p
for mod_Index in range(0, MODULUS_COUNT):
	while True:
		q_factor = int(random_prime_sympy(CIPHER_BITS))
		if q_factor not in q_base and mod_inverse(q_factor, p) == 1:
			q_base.append(q_factor)
			break

# Find new unique ciphertext moduli such that each moduli is 1 mod p
for mod_Index in range(0, SWITCH_MODULUS_COUNT):
	while True:
		q_factor = int(random_prime_sympy(SWITCH_CIPHER_BITS))
		if q_factor not in q_base and q_factor not in q_switch_base and mod_inverse(q_factor, p) == 1:
			q_switch_base.append(q_factor)
			break

q = reduce(operator.mul, q_base, 1)
q_level_down = reduce(operator.mul, q_base[1:], 1)
q_switch = reduce(operator.mul, q_switch_base, 1)


print('- plaintext modulus: ' + str(p))
print('- ciphertext modulus base: ' + str(q_base) + " = " + str(q))
print('- ciphertext level-down modulus base: ' + str(q_base) + " = " + str(q_level_down))
print('- ciphertext arbitrarily switched modulus base: ' + str(q_switch_base) + " = " + str(q_switch))
print()

A = centered_mod(np.array([random.randrange(q) for _ in range(N)], dtype=object), q)
M = centered_mod(np.array([random.randrange(p) for _ in range(N)], dtype=object), p)
E = centered_mod(np.array([random.randrange(int(10)) for _ in range(N)], dtype=object), int(10))
S = centered_mod(np.array([random.randrange(3) for _ in range(N)], dtype=object), q)

S -= 1

### Encrypt
B = poly_reduce(M + p*E - poly_reduce(np.convolve(A, S), q), q)


### Decrypt
M_decrypted = centered_mod(poly_reduce((poly_reduce(np.convolve(A, S), q) + B), q), p)


### Modulus Level Down
A_level_down = ModSwitchBGV(A, p, q, q_level_down)
B_level_down = ModSwitchBGV(B, p, q, q_level_down)
M_level_down_decrypted = centered_mod(poly_reduce(poly_reduce(np.convolve(A_level_down, S), q_level_down) + B_level_down, q_level_down), p)

### Modulus Switch
A_switch = ModSwitchBGV(A, p, q, q_switch)
B_switch = ModSwitchBGV(B, p, q, q_switch)
M_switch_decrypted = centered_mod(poly_reduce(poly_reduce(np.convolve(A_switch, S), q_switch) + B_switch, q_switch), p)



print("M_decrypted[0] = ", M_decrypted[0]) 
print("M_level_down_decrypted[0] = ", M_level_down_decrypted[0])
print("M_switch_decrypted[0] = ", M_switch_decrypted[0])
print()

print("M_decrypted == M_level_down_decrypted?", (M_level_down_decrypted == M_decrypted).all())
print("M_decrypted == M_switch_decrypted?", (M_switch_decrypted == M_decrypted).all())

