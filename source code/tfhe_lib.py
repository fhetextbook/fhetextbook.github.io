from common import round_div
import random
import numpy as np
import math

GSW_BASE = 1 << 5 # Base 32
GSW_L = 6         # to cover Q=2^30 (32^6 = 2^30

def tfhe_key_generate(K: int):
	return np.array(np.random.randint(0, 2, size=K), dtype=object)

def tfhe_encrypt(m: int, K, P, Q, s, scale):
	"""
	Returns ciphertext (vec{a}, b):
	  a <- uniform in Z_q^k
	  e <- small
	  b = -a*s + e + scale*m	(mod q)
	"""
	bound = 1
	m_scaled = scale * (m % P)
	a = np.array([random.randrange(Q) for _ in range(K)], dtype=object)
	e = random.randint(-bound, bound)
	b = (-np.dot(a, s) + m_scaled + e) % Q
	return (a, b)



def tfhe_decrypt(a, b, P, Q, s, scale, expected_dec = None):
	"""
	Decrypt:
	  v = b + a*s  (mod x^n+1, q)  ~= scale*m + e
	  m_hat = round(v/scale)
	  (optionally centered mod p if p was stored by fhe_encrypt)
	"""
	v = (np.dot(a, s) + b) % Q 
	m_hat = round_div(v, scale) % P # for the case that the result is negative 

	if expected_dec is not None:
		e = (v - expected_dec*scale) % Q
		if e > Q // 2:
			e -= Q
		max_abs = abs(e)
		e_bits = 1 if max_abs == 0 else math.floor(math.log2(abs(max_abs))) + 1
	return m_hat, e_bits


def tfhe_encrypt_without_scaling(phase_val: int, K, Q, s):
    """
    Encrypts a specific integer phase directly (no scaling, no mod P).
    Used for Circuit Bootstrapping where target phases are raw integers 
    (m * s_i * base^j).
    
    b = -a*s + phase + e
    """
    bound = 1
    a = np.array([random.randrange(Q) for _ in range(K)], dtype=object)
    e = random.randint(-bound, bound)
    
    # Encrypt the raw phase_val directly
    b = (-np.dot(a, s) + phase_val + e) % Q
    return (a, b)



def tfhe_add_cipher_cipher(a1, b1, a2, b2, Q):
	a3 = (a1 + a2) % Q
	b3 = (b1 + b2) % Q
	return (a3, b3)

def tfhe_add_cipher_plain(a1, b1, m2, Q, scale):
	b3 = (b1 + m2 * scale) % Q
	return (a1, b3) 

def tfhe_mult_cipher_plain(a1, b1, m2, Q):
	a3 = (a1 * m2) % Q
	b3 = (b1 * m2) % Q 
	return (a3, b3)


def tfhe_mult_cipher_cipher(a1, b1, a2, b2, K, P, Q, scale, relin_key):
	 """
	 TFHE (LWE) ct x ct multiplication with relinearization.
	 """
	 # 1. Compute Tensor Product (Raw)
	 b_raw = int(b1) * int(b2)
	 a_lin_raw = (int(b1) * a2 + int(b2) * a1)
	 a_quad_raw = np.outer(a1, a2).flatten()

	 # 2. Rescale
	 b_scaled = round_div(b_raw, scale) % Q
	 a_lin_scaled = np.array([round_div(x, scale) for x in a_lin_raw], dtype=object) % Q
	 a_quad_scaled = np.array([round_div(x, scale) for x in a_quad_raw], dtype=object) % Q
	 
	 # 3. Relinearize Quadratic Term
	 # The term is + <a_quad_scaled, s_old>.
	 # We want to preserve this positive relationship.
	 # Input to keyswitch should be the coefficients themselves.
	 a_in_sw = a_quad_scaled  
	 b_in_sw = 0
	 
	 a_sw, b_sw = tfhe_keyswitch(a_in_sw, b_in_sw, relin_key)
	 
	 # 4. Combine
	 # Final = (b_scaled + b_sw) + < (a_lin + a_sw), s >
	 a_final = (a_lin_scaled + a_sw) % Q
	 b_final = (b_scaled + b_sw) % Q
	 
	 return (a_final, b_final)


def unsigned_digits_len(Q: int, base: int) -> int:
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


def unsigned_gadget_decompose(x, base, L):
	"""
	Decomposes integer x into STANDARD base-B digits (Unsigned).
	x = sum(d_i * base^i) where d_i in [0, base-1]
	"""
	digits = []
	val = int(x)
	
	for _ in range(L):
		# Standard modulo gives range [0, base-1]
		d = val % base
		digits.append(d)
		
		# Standard integer division
		val = (val - d) // base
		
	return digits


def tfhe_keyswitch_keygen(s_old: np.array, s_new: np.array, K: int, Q: int, base: int = 1 << 10, bound: int = 1):
	"""
	Generates LWE Key Switching Key (KSK) from s_old (k-dim) to s_new (K-dim).
	
	Structure:
	  KSK is a list of size len(s_old).
	  Each element KSK[i] is a list of size L (levels).
	  KSK[i][j] is an LWE encryption of (s_old[i] * base^j) under s_new.
	"""
	k_old = len(s_old)
	L = unsigned_digits_len(Q, base) # from common.py
	
	ksk = []
	
	for i in range(k_old):
		ksk_i = []
		base_pow = 1
		
		for j in range(L):
			# Message to encrypt: s_old[i] * base^j
			# This is a scalar value added to the phase.
			msg = (int(s_old[i]) * base_pow) % Q
			
			# Encrypt msg under s_new
			# LWE sample: (A, b) where b = <A, s_new> + e + msg
			A = np.array([random.randrange(Q) for _ in range(K)], dtype=object)
			e = random.randint(-bound, bound)
			b = (-np.dot(A, s_new) + e + msg) % Q

			ksk_i.append((A, b))
			base_pow *= base
			
		ksk.append(ksk_i)
		
	return {"base": base, "L": L, "K": K, "Q": Q, "ksk": ksk}



def tfhe_keyswitch(a: np.array, b: int, ksk_bundle):
	"""
	Switches LWE ciphertext (a, b) from s_old to s_new.
	
	Decryption logic:
	  Phase = b - <a, s_old>
			= b - sum(a_i * s_old[i])
			= b - sum( sum(d_{i,j} * base^j) * s_old[i] )
			= b - sum( sum(d_{i,j} * (base^j * s_old[i])) )
			
	  We replace (base^j * s_old[i]) with its encryption KSK[i][j].
	"""
	base = int(ksk_bundle["base"])
	L	= int(ksk_bundle["L"])
	K	= int(ksk_bundle["K"]) # Dimension of s_new
	Q	= int(ksk_bundle["Q"])
	ksk  = ksk_bundle["ksk"]	# dim: n_old x L
	
	k_old = len(a)
	
	# Initialize result (a', b')
	# Starting with (0, b)
	a_new = np.zeros(K, dtype=object)
	b_new = b
	
	# Iterate over elements of old 'a' vector
	for i in range(k_old):
		# Decompose coefficient a[i]
		digits = unsigned_gadget_decompose(a[i], base, L)
		
		for j in range(L):
			d = digits[j]
			if d == 0:
				continue
				
			# Retrieve key switching key parts
			A_ksk, b_ksk = ksk[i][j]
			
			# Subtract d * KSK from accumulator
			# (a', b') -= d * (A_ksk, b_ksk)
			# a' -= d * A_ksk
			# b' -= d * b_ksk
			
			term_a = (d * A_ksk) % Q
			term_b = (d * b_ksk) % Q
			
			a_new = (a_new + term_a) % Q
			b_new = (b_new + term_b) % Q
			
	return (a_new, b_new)


def tfhe_circuit_bootstrap(ct_lwe, bk, K, P, Q, scale):
    """
    Converts LWE(m) -> GSW(m) using Circuit Bootstrapping.
    
    GSW(m) consists of (K+1) blocks.
    Block i corresponds to 's_i'.
    Each block contains L LWE ciphertexts.
    The (j)-th ciphertext in block i encrypts: m * s_i * base^j.
    
    For the last block (corresponding to 'b'), we encrypt: m * base^j.
    
    Note: We do not know s_i. But standard GSW construction 
    effectively requires us to generate ciphertexts that act like 
    m * Gadget_Row.
    
    In TFHE Circuit Bootstrapping:
    For each target GSW row (which corresponds to a scaling factor C),
    we run LWE Bootstrapping (PBS) on ct_lwe using a LUT that maps:
    x -> x * C.
    """
    gsw_ct = []
    
    # Dimensions of the LWE vector we are operating ON (K+1 components)
    # The GSW ciphertext must match the dimensions required by External Product.
    # If External Product takes (a_0...a_{K-1}, b), we need K+1 blocks.
    
    # We iterate 0..K. 
    # Indices 0..K-1 correspond to a_i (associated with -s_i in decryption)
    # Index K corresponds to b (associated with 1 in decryption)
    
    # Wait: The user's formula for External Product uses Decomp(ct_i).
    # To get LWE(m1*m2), if ct1 decrypts to m1 via <ct1, (-s, 1)>,
    # and GSW(m2) decrypts to m2 via <GSW, (-s, 1)>,
    # Then ExtProd decrypts to m1*m2.
    
    # Correct GSW Construction via Circuit Bootstrap:
    # We need to generate LWE encryptions of (m * H_{i,j}) where H is the gadget matrix.
    # H_{i,j} corresponds to the value bit `base^j` at position `i`.
    # Ideally, H_{i,j} is related to the secret key s_i, but in Circuit Bootstrap
    # we generate encryptions of 'm * base^j' and rely on the external product
    # structure to route them correctly.
    
    # Actually, for the User's specific "Cipher-Cipher" multiplication:
    # We just need to generate LWE encryptions of (m * base^j) for the last column
    # (the 'b' part) and 0 for the 'a' parts? 
    # No, Circuit Bootstrap generates a full valid GSW ciphertext.
    # It requires bootstrapping keys that encrypt the bits of the secret key `s`.
    # Since we cannot implement the full RGSW/GLWE keygen here (too large),
    # we will implement a FUNCTIONAL SIMULATION of Circuit Bootstrap 
    # that computes the correct LWEs directly if 's' is provided (for demo),
    # OR we implement the programmable bootstrap logic.
    
    # FOR THIS DEMO: We will implement the PBS-based generation logic.
    # We assume 'bk' contains the necessary material to compute PBS.
    
    base = GSW_BASE
    L = GSW_L
    
    for i in range(K + 1):
        block = []
        for j in range(L):
            scaling_factor = (base**j)
            
            # If i < K, this corresponds to a_i part.
            # In standard GSW, row i encrypts s_i * m * base^j.
            # We need a PBS that outputs LWE(s_i * m * base^j).
            # This requires a "Bootstrapping Key" that can select on s_i.
            # Since standard PBS only selects on 'm' (phase of ct_lwe), 
            # we need specific keys.
            
            # SIMPLIFICATION FOR DEMO:
            # We will perform the mathematical equivalent operation directly
            # to show the "Circuit Bootstrap" flow without the massive overhead
            # of GLWE/FFT implementation in Python.
            
            # Ideally: lwe_out = tfhe_bootstrap(ct_lwe, lut_func=lambda x: x * scaling_factor * s_i)
            # But we don't know s_i inside this function unless passed.
            
            # We will use a helper that simulates the PBS output for the sake of the "Cipher-Cipher" demo.
            # In a real C++ lib, this uses the TGsw samples in the Bootstrapping Key.
            
            if i < K:  # The first k blocks in GSW
                # Need LWE( m * s_i * base^j )
                # Real Circuit Bootstrap uses a 'Select' Key (LWE encryptions of s_i).
                # We will simulate this by passing 's' in the 'bk' for this Python prototype.
                s_val = bk['s'][i] # Cheating slightly for prototype size limits
                
                # LUT: y = x * s_val * scaling_factor
                # We assume LWE encrypts 'm' scaled by 'scale'.
                # We want LWE of result scaled by 'scale'.
                
                def lut_func(x):
                    # x is the message m
                    val = (x * s_val * scaling_factor)
                    return val
                
                # Run PBS (Simulated or Real)
                lwe_sample = tfhe_mock_pbs(ct_lwe, lut_func, K, P, Q, scale, bk['s'])
                
            else:   # the last block in GSW
                # i == K (The constant term 1)
                # Need LWE( m * 1 * base^j )
                def lut_func(x):
                    val = (x * scaling_factor)
                    return val
                
                lwe_sample = tfhe_mock_pbs(ct_lwe, lut_func, K, P, Q, scale, bk['s'])
            
            block.append(lwe_sample)
        gsw_ct.append(block)
        
    return gsw_ct



def tfhe_mock_pbs(ct_in, lut_func, K, P, Q, scale, s):
    """
    Performs a Programmable Bootstrapping (PBS) functionally.
    In a real library, this uses Blind Rotate + Sample Extract.
    Here, we decrypt, apply LUT, and re-encrypt to allow the
    higher-level 'Cipher-Cipher' logic to be tested.
    """
    # 1. Decrypt to get the underlying message m
    a, b = ct_in
    phase = (b + np.dot(a, s)) % Q
    m_noisy = round_div(phase, scale) % P
    
    # 2. Apply LUT to get the new TARGET PHASE
    # Note: We do NOT mod P here. The result is a raw phase in Z_Q.
    target_phase = lut_func(m_noisy) % Q
    
    # 3. Re-encrypt the raw phase directly (Bypassing tfhe_encrypt scaling)
    return tfhe_encrypt_without_scaling(target_phase, K, Q, s)


# ==========================================
#        Main Multiplication Wrapper
# ==========================================

def tfhe_mult_cipher_cipher_circuit_bootstrapping(a1, b1, a2, b2, K, P, Q, scale, bk):
    """
    Cipher-Cipher Multiplication using Circuit Bootstrapping.
    ct1: LWE(m1)
    ct2: LWE(m2)
    Returns: LWE(m1 * m2)
    
    Method:
    1. Convert ct2 -> GSW(m2) using Circuit Bootstrap.
    2. Compute ct1 * GSW(m2) using External Product.
    """

    L = unsigned_digits_len(Q, GSW_BASE)

    ct1 = (a1, b1)
    ct2 = (a2, b2)

    # 1. Circuit Bootstrap: LWE(m2) -> GSW(m2)
    # Note: We pass 'bk' which contains the secret key 's' for our mock PBS.
    gsw_m2 = tfhe_circuit_bootstrap(ct2, bk, K, P, Q, scale)
    
    # 2. External Product: LWE(m1) * GSW(m2) -> LWE(m1 * m2)
    # Formula: sum( <Decomp(ct1_i), GSW_row_i> )
    res_ct = tfhe_external_product(ct1, gsw_m2, K, Q, base=GSW_BASE, L=GSW_L)
    
    return res_ct



def tfhe_external_product(ct_lwe, ct_gsw, K, Q, base=GSW_BASE, L=GSW_L):
    """
    Computes LWE * GSW -> LWE.
    Formula: sum( < Decomp(ct_i), GSW_i > )
    """
    a, b = ct_lwe
    # Combine a and b into a single vector for iteration: (a_0, ... a_{k-1}, b)
    # This matches the user notation: ct = (ct_0, ... ct_k)
    full_ct = np.concatenate([a, [b]])
    
    # Initialize result accumulator (0, 0)
    res_a = np.zeros(K, dtype=object)
    res_b = 0
    
    # Iterate through each component of the LWE ciphertext
    # ct_gsw is a list of 'Lev' blocks. ct_gsw[i] corresponds to full_ct[i].
    for i in range(len(full_ct)):
        # Decompose the scalar component full_ct[i]
        digits = unsigned_gadget_decompose(full_ct[i], base, L)
        
        # Access the i-th block of the GSW ciphertext
        # This block is a list of L ciphertexts (Lev_0, ... Lev_{L-1})
        lev_block = ct_gsw[i]
        
        for j in range(L):
            d = digits[j]
            if d == 0: continue
            
            # The j-th element of the Lev block is an LWE ciphertext
            lev_a, lev_b = lev_block[j]
            
            # Add d * Lev_j to result
            term_a = (d * lev_a) % Q
            term_b = (d * lev_b) % Q
            
            res_a = (res_a + term_a) % Q
            res_b = (res_b + term_b) % Q
            
    return (res_a, res_b)



def tfhe_bootstrap_keygen(s: np.array):
    """
    Generates the Bootstrapping Key (BK) required for Circuit Bootstrapping.
    
    In a FULL TFHE implementation, this would generate RGSW encryptions 
    of the secret key bits (for Blind Rotation).
    
    In this FUNCTIONAL PROTOTYPE (Mock PBS), this packages the secret key 's'
    to allow the tfhe_mock_pbs function to simulate the bootstrapping noise 
    refresh and LUT application.
    """
    bk = {
        's': s,
        'type': 'mock_bk',
        'info': 'Contains secret key for functional PBS simulation'
    }
    return bk



