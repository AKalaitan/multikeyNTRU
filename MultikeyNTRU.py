from NTRU_lib import *
import time

#module for testing and benchmarking the scheme
#K - is the power for n; B - is the bound; N - number of the parties; M - array of the messages 
def test(k,B,N,m):
	n = 2**k
	sec = int(math.sqrt(n))
	print "%d-bit security level" %sec
	bit_len = (2*n*B*B*(2*n*B+1)*(2*B+1))**N
	
	#generate a prime number q, for Rq
	q = int(number.getPrime(bit_len.bit_length()+1))

	#G - is the Ideal
	G = []
	#f - is the private key
	f = []
	#h - is the public key
	h = []
	# s and e - are 'noise' polynomials
	s = []
	e = []

	###########################################
	#Generate N keys:
	start_time = time.time()
	#both Gaussian and seed_generated distributions work, specify which one you need
	#for the Gaussian distribution the scipy package required
	###########################################
	for i in range(0,N):
		#G1,f1,h1,s1,e1 = keyGen(n,B,q) # take polynomials from Gaussian distribution
		G1,f1,h1,s1,e1 = keyGen_from_random_distr(n,B,q) # take B-bounded polynomials with random coefficients
		G.append(G1)
		f.append(f1)
		h.append(h1)
		s.append(s1)
		e.append(e1)

	print "Secret key - bit-size = %d" %(bit_size(f[0])) 
	print "Public key - bit-size = %d" %(bit_size(h[0])) 

	print "Key gen time :"
	print("--- %s seconds ---" % (time.time() - start_time))

	start_time = time.time()

	###########################################
	# Encrypt N messages:
	start_time = time.time()
	###########################################
	c =[]
	
	for i in range(0,N):
		c.append(enc(m[i],h[i],s[i],e[i],q)) 

	print "Ciphertext - bit-size = %d" %(bit_size(c[0])) 

	print "Encryption time :"
	print("--- %s seconds ---" % (time.time() - start_time))

	
	#############################################
	#Homomorphical evaluation of the Summation operation
	start_time = time.time()
	#############################################

	c_add = [0]
	for i in range(0,N):
		c_add = poly_add(c_add,c[i],q)

	print "Joint ciphertext of addition - bit-size = %d" %(bit_size(c_add)) 

	print "Addition of %d ciphertexts time :" %N
	print("--- %s seconds ---" % (time.time() - start_time))

	#############################################
	#Homomorphical evaluation of the Multiplication operation
	start_time = time.time()
	#############################################
	
	c_mult = [1]
	for i in range(0,N):
		c_mult = poly_mult(c_mult,c[i],q)

	print "Joint ciphertext of multiplication - bit-size = %d" %(bit_size(c_mult)) 

	print "Multiplication of %d ciphertexts time :" %N
	print("--- %s seconds ---" % (time.time() - start_time))

	#############################################
	#Combining of the keys 
	#This step is made for simplicity, 
	#but we still can consequently multiply the joint ciphertext by the secret keys, with the correct result at the end
	start_time = time.time()
	#############################################

	f_joint = [1]
	for i in range(0,N):
		f_joint = poly_mult(f_joint,f[i],q)
	
	print "Combined secret key - bit-size = %d" %(bit_size(f_joint)) 

	print "Combined key computation time :"
	print("--- %s seconds ---" % (time.time() - start_time))

	

	#############################################
	#Decryption of the C_add
	start_time = time.time()
	#############################################

	check_add =0
	for i in range(0,N):
		print "message[%d] = %d" %(i+1,m[i][0])
		check_add ^=m[i][0]

	for i in range(0,N):
		print m[i], " XOR ",
	print " = ", check_add


	print "Decryption of C_add:"
	print dec(c_add,f_joint,G[0],q)

	print "Decryption of C_add time :"
	print("--- %s seconds ---" % (time.time() - start_time))

	#############################################
	#Decryption of the C_mult
	start_time = time.time()
	#############################################

	check_mult = 1
	for i in range(0,N):
		print "message[%d] = %d" %(i+1,m[i][0])
		check_mult *=m[i][0]

	for i in range(0,N):
		print m[i], " x ",
	print " = ", check_mult


	print "Decryption of C_mult:"
	print dec(c_mult,f_joint,G[0],q)

	print "Decryption of C_mult time :"
	print("--- %s seconds ---" % (time.time() - start_time))

	
#############################################
#Main function# Put any values you want below
#############################################
k = 5  # n = 2**k, so n is the power of 2; n - is the maximum degree of the polynomials over the ring Rq
B = 6  # B-bound for the distribution; all the coefficients are less than B-bound
N = 7  # The number of parties
m = [[1],[0],[1],[1],[0],[1],[1]] # The messages of each party
test(k,B,N,m)			  # Run encryption

