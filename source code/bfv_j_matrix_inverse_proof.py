import numpy as np
from sympy import nextprime, mod_inverse, primitive_root
import copy

M = 32

N = int(M / 2)
P = nextprime(M)
while (P - 1) % M != 0:
   P = nextprime(P)

mult_group_generator = primitive_root(P)  # prim_root is a number such that Order(prim_root) = totient(p) mod p
primitive_root = pow(mult_group_generator, (P - 1) // M, P)


print("N: " + str(N))
print("P: " + str(P))
print("Primitive Root of Unity: " + str(primitive_root))
print()

for i in range(int(M/4)):
	print("J(" + str(i) + ") = 5^{" + str(i) + "} mod 2*" + str(N) + " = " + str((5 ** (i)) % M) )
print()

for i in reversed(range(int(M/4))):
	print("J_*(" + str(i) + ") = -5^{" + str(i) + "} mod 2*" + str(N) + " = " + str((-5 ** (i)) % M) )
print()
print()
WT = np.zeros((N, N), dtype = object)

print('---------------------')

for i in range(N):
	for j in range(N):
		if i < N / 2:
			WT[i][j] = (((5 ** i) % M) * j ) % M
		else:
			WT[i][j] = (((-5 ** (N - i - 1)) % M) * j ) % M
			#WT[i][j] = (((-5 ** i) % M) * j ) % M
	print()
W = WT.T

print("<W Matrix>")
print(W)
print("<W^T Matrix>")
print(WT)

W_modified = copy.deepcopy(W)
WT_modified = copy.deepcopy(WT)

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

print("<W-modified Matrix>")
print(W_modified)
print("<WT-modified Matrix>")
print(WT_modified)

Sum_Matrix = np.zeros((N, N), dtype = object)

for i in range(N):
	for j in range(N):
		print("<WT*W(" + str(i) + "," + str(j) + ") Exponent>")
		arr = []
		summation = 0
		for k in range(N):
			arr.append((WT[i][k] + W[k][j]) % M)
			summation += primitive_root ** (arr[len(arr) - 1])
		print("Exponent List: " + str(arr))
		print("Summation: " + str(summation % P)) 
		print()
		Sum_Matrix[i][j] = summation % P

print("<W^T*W Actual Summation>")
print(Sum_Matrix)

Modified_Sum_Matrix = np.zeros((N, N), dtype = object)

for i in range(N):
	for j in range(N):
		print("<Modified WT*W(" + str(i) + "," + str(j) + ") Exponent>")
		arr = []
		summation = 0
		for k in range(N):
			arr.append((WT_modified[i][k] + W_modified[k][j]) % M)
			summation += primitive_root ** (arr[len(arr) - 1])
		print("Exponent List: " + str(arr))
		print("Summation: " + str(summation % P)) 
		print()
		Modified_Sum_Matrix[i][j] = summation % P

print("<Modified W^T*W Actual Summation>")
print(Modified_Sum_Matrix)
#print(vector**2)
#print(added)
