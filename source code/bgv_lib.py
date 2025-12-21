from common import *


def fhe_encrypt(m: Polynomial, N, P, Q_BASE, s, scale):
	"""
	Returns ciphertext (a, b):
	 a <- uniform in R_q
	 e <- small
	 b = -a*s + scale*e + m	(mod x^n+1, q)
	"""
	Q = math.prod(Q_BASE)

	m_mod = poly_mod(m, N, P)  # coeffs in [0,p)
	a = sample_uniform_poly(N, Q)
	e = sample_small_poly(N, bound=1)
	e_scaled = Polynomial(np.array([c * scale for c in e.coef], dtype=object))

	b = poly_add(poly_add(-poly_mult(a, s, N, Q), e_scaled, N, Q),
				m_mod, N, Q)

	a = poly_mod(a, N, Q)
	b = poly_mod(b, N, Q)
	return (a, b)



def fhe_decrypt(a: Polynomial, b: Polynomial, N, P, Q_BASE, s, scale, expected_dec: Polynomial = None) -> Polynomial:
	"""
	Decrypt:
	 v = b + a*s  (mod x^n+1, q)  ~= m + e*scale
	 m_hat = round(v/scale)
	 (optionally centered mod p if p was stored by fhe_encrypt)
	"""
	Q = math.prod(Q_BASE)

	v = poly_add(b, poly_mult(a, s, N, Q), N, Q)
	m_hat = poly_mod(v, N, P)
	print(m_hat)
	if expected_dec is not None:
		noise_poly = poly_sub(v, expected_dec, N, Q)
		e = Polynomial(np.array([round_div(c, P) for c in noise_poly.coef], dtype=object))
		max_abs = float(np.max(np.abs(e.coef))) if len(e.coef) else 0.0
		e_bits = 1 if max_abs == 0 else math.floor(math.log2(abs(max_abs))) + 1
	return m_hat, e_bits


def bgv_mod_switch_poly(poly, N, P, Q_cur, Q_new):
	"""
	Performs BGV Modulus Switching from Q_cur to Q_new.
	Ensures c' = c mod P.
	Logic ported from bgv.py's ModSwitchBGV using integer arithmetic.
	"""
	# 1. Prepare arrays (ensure dtype=object for large integers)
	coeffs = np.array(poly.coef, dtype=object)
	
	# 2. Scale: c_scaled = round(c * Q_new / Q_cur)
	#	We implement vectorized round_div logic here:
	#	round(x/d) = (x + d//2) // d		 if x >= 0
	#					= -((-x + d//2) // d)	if x < 0
	
	num = coeffs * Q_new
	half_Q = Q_cur // 2
	
	# Use np.where to handle positive and negative coefficients efficiently
	c_scaled = np.where(
		 num >= 0,
		 (num + half_Q) // Q_cur,
		 -((-num + half_Q) // Q_cur)
	)
	
	# 3. Compute error term (epsilon)
	#	epsilon = c * Q_new - c_scaled * Q_cur
	epsilon = num - (c_scaled * Q_cur)
	
	# 4. Compute Correction H
	#	H = (epsilon * inv(Q_cur)) mod P
	inv_Q_cur = mod_inv(Q_cur, P)
	if inv_Q_cur == -1:
			raise ValueError(f"Q_cur={Q_cur} is not coprime to P={P}")
			
	h = (epsilon * inv_Q_cur) % P
	h = centered_mod_arr(h, P) # Vectorized centered mod from common.py
	
	# 5. Final Result: scaled + H (mod Q_new)
	result = c_scaled + h
	new_coeffs = centered_mod_arr(result, Q_new)
		 
	return Polynomial(np.array(new_coeffs, dtype=object))


def fhe_add_cipher_cipher(a1, b1, a2, b2, N, Q_BASE):
	Q = math.prod(Q_BASE)
	a3 = poly_add(a1, a2, N, Q)
	b3 = poly_add(b1, b2, N, Q)
	return (a3, b3)



def fhe_add_cipher_plain(a1, b1, p2, N, Q_BASE, scale):
	Q = math.prod(Q_BASE)
	b3 = poly_add(b1, p2, N, Q)
	return (a1, b3) 

def fhe_mult_cipher_plain(a1, b1, p2, N, Q):
	a3 = poly_mult(a1, p2, N, math.prod(Q))
	b3 = poly_mult(b1, p2, N, math.prod(Q))
	return (a3, b3)


def fhe_rotate(a1, b1, N, Q_BASE, rotation_offset, galois_rotation_key):
	Q = math.prod(Q_BASE)
	# 1) Automorphism on ciphertext components
	a_rot = poly_mod(poly_rotate(a1, Q, N, rotation_offset), N, Q)
	b_rot = poly_mod(poly_rotate(b1, Q, N, rotation_offset), N, Q)

	# 2) Key-switch from s_rot -> s (cache by (offset,N,Q,s)) by gadget decomposition
	a_out, b_out = fhe_keyswitch(a_rot, b_rot, galois_rotation_key)
	return (a_out, b_out)




def fhe_mult_cipher_cipher(a1, b1, a2, b2, N, P, Q_BASE, scale, relin_key):
	"""
	BGV ct×ct with gadget decomposition relinearization.
	"""
	# Use the plaintext modulus from the relin key bundle if you store it elsewhere.
	Q = math.prod(Q_BASE)

	# Ensure inputs are centered in R_Q
	a1 = poly_mod(a1, N, Q); b1 = poly_mod(b1, N, Q)
	a2 = poly_mod(a2, N, Q); b2 = poly_mod(b2, N, Q)

	# ---- Raw products in big integers (no mod Q yet) ----
	d0_raw = poly_mult(b1, b2, N, Q)  # b1*b2
	d1_raw = poly_add(
		poly_mult(b1, a2, N, Q),
		poly_mult(a1, b2, N, Q),
		N, Q
	)  # b1*a2 + a1*b2
	d2_raw = poly_mult(a1, a2, N, Q)  # a1*a2

	# ---- BGV rescale: round(P/Q * ·) and reduce mod Q ----
	q_last = Q_BASE.pop()
	Q_new = Q // q_last

	d0 = bgv_mod_switch_poly(d0_raw, N, P, Q, Q_new)
	d1 = bgv_mod_switch_poly(d1_raw, N, P, Q, Q_new)
	d2 = bgv_mod_switch_poly(d2_raw, N, P, Q, Q_new)


	# ---- Relinearize the s^2 term using the gadget-decomp keyswitch ----
	# Under s^2, (a_old=d2, b_old=0) has phase = d2*s^2
	a_old = d2
	b_old = poly_zero(N)
	a_sw, b_sw = fhe_keyswitch(a_old, b_old, relin_key)

	# Combine into 2-term ciphertext under s
	a3 = poly_add(d1, a_sw, N, Q_new)
	b3 = poly_add(d0, b_sw, N, Q_new)
	return (a3, b3)



def fhe_bootstrapping(a1, b1, N, Q):
	return 0
