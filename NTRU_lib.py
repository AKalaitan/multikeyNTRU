from poly_lib import *
from Crypto.Util import number

#The last step of decryption
def mod_2(f):
	for i in range(0,len(f)):
		if f[i]%2 != 0 : return False
	return True

#Generate the keys
def keyGen(n,B_bound,q):
	G = [1]
	for i in range(0,n-1):
		G.append(0)
	G.append(1)
	n = n/4 # while the degree of f_inverse is n, we do not need to use "noise" polynomials with the same degree

	finv = [0]

	while finv == [0] or len(finv) == 1:
		f = take_f_from_gaussian_distribution(n,B_bound,q) #take_f_from_cyclo_distribution(n,q)
		g = gaussian(n,B_bound)
		s = gaussian(n,B_bound)
		e = gaussian(n,B_bound)

		while (check_quo(poly_mult(g,s,q),q) == False) :
			g = gaussian(n,B_bound)
			s = gaussian(n,B_bound)

		while (check_quo(e,q) == False) :
			e = gaussian(n,B_bound)

		gcd,finv,ginv = xgcd(f,G,q)
	
	h = poly_mult([2],poly_mult(finv,g,q),q)
	e = poly_mult([2],e,q)
	return G,f,h,s,e

def keyGen_from_random_distr(n,B_bound,q):
	G = [1]
	for i in range(0,n-1):
		G.append(0)
	G.append(1)
	n = n/4 # while the degree of f_inverse is n, we do not need to use "noise" polynomials with the same degree

	finv = [0]

	while finv == [0] or len(finv) == 1:
		f = take_f_from_seed_gen_distribution(n,B_bound,q) #take_f_from_cyclo_distribution(n,q)
		g = randomKeyPoly(n,B_bound)
		s = randomKeyPoly(n,B_bound)
		e = randomKeyPoly(n,B_bound)

		while (check_quo(poly_mult(g,s,q),q) == False) :
			g = randomKeyPoly(n,B_bound)
			s = randomKeyPoly(n,B_bound)

		while (check_quo(e,q) == False) :
			e = randomKeyPoly(n,B_bound)

		gcd,finv,ginv = xgcd(f,G,q)
	
	h = poly_mult([2],poly_mult(finv,g,q),q)
	e = poly_mult([2],e,q)
	return G,f,h,s,e


#Encryption
def enc(m,h,s,e,q):
	return poly_add(poly_add(poly_mult(h,s,q),e,q),m,q)

#Decryption
def dec(c,f,G,q):
	m_dec = poly_mod(poly_mult(c,f,q),G,q)
	if mod_2(m_dec) == True : return 0
	else : return 1

#Returns the bit-size of an object
def bit_size(f):
	count = 0
	for i in range(0,len(f)):
		count +=f[i].bit_length()
	return count
