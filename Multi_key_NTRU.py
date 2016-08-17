from poly_lib import *
from Crypto.Util import number
import time

def mod_2(f):
	for i in range(0,len(f)):
		if f[i]%2 != 0 : return False
	return True

def keyGen(n,B_bound,q):
	G = [1]
	for i in range(0,n-1):
		G.append(0)
	G.append(1)
	n = n/4 #check it! 

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

def enc(m,h,s,e,q):
	return poly_add(poly_add(poly_mult(h,s,q),e,q),m,q)

def dec(c,f,G,q):
	m_dec = poly_mod(poly_mult(c,f,q),G,q)
	if mod_2(m_dec) == True : return 0
	else : return 1

def bit_size(f):
	count = 0
	for i in range(0,len(f)):
		count +=f[i].bit_length()
	return count
#k = 5 # 5 tested
k = 9 #22-bit secure parameter
n = 2**k
B_bound = 6
B = B_bound

B = B_bound
bit_len = 2*n*B*B*(2*n*B+1)*(2*B+1)
q = int(number.getPrime(bit_len.bit_length()+1))

print "22-bit security level are in testing"
print "q - bit-size = %d" %(q.bit_length())

#key 1
###########################################
start_time = time.time()
G,f1,h1,s1,e1 = keyGen(n,B_bound,q)

print "Secret key 1 - bit-size = %d" %(bit_size(f1)) 
print "Public key 1 - bit-size = %d" %(bit_size(h1)) 


print "Key gen 1 time :"
print("--- %s seconds ---" % (time.time() - start_time))


m0_1 = [0]

start_time = time.time()

c0_1 = enc(m0_1,h1,s1,e1,q)

print "Ciphertext 1 - bit-size = %d" %(bit_size(c0_1)) 


print "Encryption 1 time :"
print("--- %s seconds ---" % (time.time() - start_time))

#print c0_1
#print dec(c0_1,f1,G,q)
#key 2
############################################
start_time = time.time()

G,f2,h2,s2,e2 = keyGen(n,B_bound,q)

print "Secret key 1 - bit-size = %d" %(bit_size(f2)) 
print "Public key 1 - bit-size = %d" %(bit_size(h2)) 

print "Key gen 2 time :"
print("--- %s seconds ---" % (time.time() - start_time))

m0_2 = [0]

start_time = time.time()

c0_2 = enc(m0_2,h2,s2,e2,q)

print "Ciphertext 2 - bit-size = %d" %(bit_size(c0_2)) 

print "Encryption 2 time :"
print("--- %s seconds ---" % (time.time() - start_time))
#print c0_2
#print dec(c0_2,f2,G,q)


#key 3
#############################################

start_time = time.time()

G,f3,h3,s3,e3 = keyGen(n,B_bound,q)

print "Secret key 1 - bit-size = %d" %(bit_size(f3)) 
print "Public key 1 - bit-size = %d" %(bit_size(h3)) 

print "Key gen 3 time :"
print("--- %s seconds ---" % (time.time() - start_time))

m0_3 = [0]

start_time = time.time()

c0_3 = enc(m0_3,h3,s3,e3,q)

print "Ciphertext 3 - bit-size = %d" %(bit_size(c0_3)) 

print "Encryption 3 time :"
print("--- %s seconds ---" % (time.time() - start_time))

#print c0_2
#print dec(c0_3,f3,G,q)

#homomorphical evaluation
#############################################
start_time = time.time()

c_add = poly_add(c0_1,c0_2,q)

#print "IN ZZ"
#print poly_mult_in_ZZ(poly_mult_in_ZZ(f1,f2),c_add)

#print "degree of F1*F2 = %d" %(len(f12))

#print "add of c1 & c2 :"
#print dec(c_add,f12,G,q)
#print "mult of c1 & c2 :"
#print dec(c_mult,f12,G,q)


#print "degree of F1*F2*F3 = %d" %(len(f123))
c_add_3 = poly_add(c_add,c0_3,q)

print "Ciphertext of addition of 3 - bit-size = %d" %(bit_size(c_add_3)) 


print "Addition of 3 time :"
print("--- %s seconds ---" % (time.time() - start_time))

start_time = time.time()

c_mult = poly_mult(c0_1,c0_2,q)
c_mult_3 = poly_mult(c_mult,c0_3,q)

print "Ciphertext of multiplication of 3 - bit-size = %d" %(bit_size(c_mult_3)) 

print "Multiplication of 3 time :"
print("--- %s seconds ---" % (time.time() - start_time))

start_time = time.time()

f12 = poly_mult(f1,f2,q)
f123 = poly_mult(f12,f3,q)

print "Combined secret key - bit-size = %d" %(bit_size(f123)) 

print "Combined key computation time :"
print("--- %s seconds ---" % (time.time() - start_time))

start_time = time.time()

print "add of c1 & c2 & c3:"
print dec(c_add_3,f123,G,q)

print "Decryption of addition of 3 time :"
print("--- %s seconds ---" % (time.time() - start_time))

start_time = time.time()

print "mult of c1 & c2 & c3:"
print dec(c_mult_3,f123,G,q)

print "Decryption of multiplication of 3 time :"
print("--- %s seconds ---" % (time.time() - start_time))




