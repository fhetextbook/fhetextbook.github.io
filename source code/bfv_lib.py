from common import *



def fhe_encrypt(m: Polynomial, N, P, Q, s, scale):
	"""
	Returns ciphertext (a, b):
	  a <- uniform in R_q
	  e <- small
	  b = -a*s + e + scale*m   (mod x^n+1, q)
	"""
	m_red = poly_mod(m, N, P)  # coeffs in [0,p)
	m_scaled = Polynomial(np.array([c * scale for c in m_red.coef], dtype=object))

	a = sample_uniform_poly(N, Q)
	e = sample_small_poly(N, bound=1)

	b = poly_add(poly_add(-poly_mult(a, s, N, Q), e, N, Q),
				 m_scaled, N, Q)

	a = poly_mod(a, N, Q)
	b = poly_mod(b, N, Q)
	return (a, b)



def fhe_decrypt(a: Polynomial, b: Polynomial, N, P, Q, s, scale, expected_dec: Polynomial = None) -> Polynomial:
	"""
	Decrypt:
	  v = b + a*s  (mod x^n+1, q)  ~= scale*m + e
	  m_hat = round(v/scale)
	  (optionally centered mod p if p was stored by fhe_encrypt)
	"""

	v = poly_add(b, poly_mult(a, s, N, Q), N, Q)
	m_hat = [round_div(c, scale) for c in v.coef]
	m_hat = [centered_mod(c, P) for c in m_hat]

	if expected_dec is not None:
		e = poly_sub(v, poly_scalar_mult(expected_dec, scale, N, Q), N, Q)
		max_abs = float(np.max(np.abs(e.coef))) if len(e.coef) else 0.0
		e_bits = 1 if max_abs == 0 else math.floor(math.log2(abs(max_abs))) + 1
	return Polynomial(np.array(m_hat, dtype=object)), e_bits

def fhe_add_cipher_cipher(a1, b1, a2, b2, N, Q):
	a3 = poly_add(a1, a2, N, Q)
	b3 = poly_add(b1, b2, N, Q)
	return (a3, b3)

def fhe_add_cipher_plain(a1, b1, p2, N, Q, scale):
	b3 = poly_add(b1, poly_scalar_mult(p2, scale, N, Q), N, Q)
	return (a1, b3) 

def fhe_mult_cipher_plain(a1, b1, p2, N, Q):
	a3 = poly_mult(a1, p2, N, Q)
	b3 = poly_mult(b1, p2, N, Q) 
	return (a3, b3)

def fhe_rotate(a1, b1, N, Q, rotation_offset, galois_rotation_key):
    # 1) Automorphism on ciphertext components
    a_rot = poly_mod(poly_rotate(a1, Q, N, rotation_offset), N, Q)
    b_rot = poly_mod(poly_rotate(b1, Q, N, rotation_offset), N, Q)

    # 2) Key-switch from s_rot -> s (cache by (offset,N,Q,s)) by gadget decomposition
    a_out, b_out = fhe_keyswitch(a_rot, b_rot, galois_rotation_key)
    return (a_out, b_out)






def fhe_mult_cipher_cipher(a1, b1, a2, b2, N, P, Q, scale, relin_key):
    """
    BFV ct×ct with gadget decomposition relinearization.
    IMPORTANT: uses BFV rescale round(P/Q * ·), NOT division by 'scale'.
    """
    # Use the plaintext modulus from the relin key bundle if you store it elsewhere.

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

    # ---- Relinearize the s^2 term using the gadget-decomp keyswitch ----
    # Under s^2, (a_old=d2, b_old=0) has phase = d2*s^2
    a_old = poly_mod(d2, N, Q)
    b_old = poly_zero(N)
    a_sw, b_sw = fhe_keyswitch(a_old, b_old, relin_key)

    # Combine into 2-term ciphertext under s
    a3 = poly_add(d1, a_sw, N, Q)
    b3 = poly_add(d0, b_sw, N, Q)
    return (a3, b3)


def fhe_bootstrapping(a1, b1, N, P, Q):
	return 0
