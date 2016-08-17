# math and scipy are used for discrete gaussian distribution 
from scipy.special import kn
import random
import math

#return gcd,f_inv,g_inv
def xgcd(f,g,q):
	f_store = []
	f_store.append(f)
	inverse = []
	while len(g) != 1:
		g, f = poly_div(f,g,q)[1], g
		f_store.append(f)
	if poly_mult([modulo_reduction(int_xgcd(g[0],q)[0],q)],g,q) == [1] : 
		gcd_poly = [1]
		inverse.append([modulo_reduction(int_xgcd(g[0],q)[0],q)])
	else: 
		return poly_xgcd(f_store[0],f_store[1],q)

	for i in range(0,len(f_store)-1):
		inverse.append(poly_div(poly_sub([1],poly_mult(f_store[len(f_store)-2-i],inverse[i],q),q),f_store[len(f_store)-1-i],q)[0])
	
	return [gcd_poly,inverse[len(inverse)-2],inverse[len(inverse)-1]]

#if gcd!=1 then use this standard algo: return c = gcd and a = inverse of f  and b = inverse of g  
def poly_xgcd(f,g,q):
	if len(g) == 1 and g == [0] :
		return [f,[1],[0]]
	else:
		c,a,b = poly_xgcd(g,poly_div(f,g,q)[1],q)
		return [c, b, poly_sub(a,poly_mult(poly_div(f,g,q)[0],b,q),q)]

# xgcd for integers
def int_xgcd(a,b):
	if b == 0:
		return [1,0,a]
	else:
		x,y,d = int_xgcd(b, a%b)
		return [y, x - (a//b)*y, d]

#divide f by g in Rq
def poly_div(f,g,q):
	if len(g) == 1:
		quo = []
		for i in range(0,len(f)):
			quo.append(0)
		for i in range(0,len(f)):
			if g[0]*(f[i]/g[0]) != f[i] : 
				temp_g = modulo_reduction((int_xgcd(g[0],q)[0]),q)
				quo[i] += (modulo_reduction((f[i]*temp_g),q))
			else: 
				quo[i] +=(modulo_reduction((f[i]/g[0]),q))
		return (quo,[0])
	
	if len(f)<len(g) : 
		#f,g = swap(f,g)
		return ([0],f)
	f_original = f[:]
	quo = []
	for i in range(0,len(f)):
		quo.append(0)
	j = 0
	while len(f)>=len(g):
		if g[j]*(f[j]/g[j]) != f[j] : 
			temp_g = modulo_reduction((int_xgcd(g[j],q)[0]),q)
			quo[len(f)-len(g)] += (modulo_reduction((f[j]*temp_g),q))
		else: 
			quo[len(f)-len(g)] +=(modulo_reduction((f[j]/g[j]),q))
		quo_step = quo[:]
		quo_step.reverse()
		quo_step = remove_lead_zero(quo_step)
		f = poly_sub(f_original , poly_mult(g,quo_step,q) ,q)
	rem = f
	quo.reverse()
	quo = remove_lead_zero(quo)
	return (quo,rem)

# compute f mod g in q
def poly_mod(f,g,q):
	rem = f[:]
	while len(rem)>=len(f):
		rem = poly_div(f,g,q)[1]
		if rem == f : return f
	return rem

# multiply f by g in Rq
def poly_mult(f,g,q):
	size = len(f)+len(g)-1
	ans = [0]
	for i in range(0,len(f)):
		size -= 1
		step_result = []
		for j in range(0,len(g)):
			step_result.append(modulo_reduction((f[i]*g[j]),q))
		while j<size:
			step_result.append(0)
			j+=1
		ans = poly_add(step_result,ans,q)
	return ans

# multiply f by g in ZZ['x']
def poly_mult_in_ZZ(f,g):
	size = len(f)+len(g)-1
	ans = [0]
	for i in range(0,len(f)):
		size -= 1
		step_result = []
		for j in range(0,len(g)):
			step_result.append(f[i]*g[j])
		while j<size:
			step_result.append(0)
			j+=1
		ans = poly_add_in_ZZ(step_result,ans)
	return ans

# add f to g in ZZ['x'] 
def poly_add_in_ZZ(f1,g1):
	f = f1[:]
	g = g1[:]
	f.reverse()
	g.reverse()
	if len(f)>len(g):
		for i in range(0,len(g)):
			f[i] += g[i]
		f.reverse()
		g.reverse()
		return f
	else:
		for i in range(0,len(f)):
			g[i] += f[i]
		f.reverse()
		g.reverse()
		return g

# add f to g in Rq 
def poly_add(f1,g1,q):
	f = f1[:]
	g = g1[:]
	f.reverse()
	g.reverse()
	if len(f)>len(g):
		for i in range(0,len(g)):
			f[i] += g[i]
			f[i] = modulo_reduction(f[i],q)
		f.reverse()
		g.reverse()
		return f
	else:
		for i in range(0,len(f)):
			g[i] += f[i]
			g[i] = modulo_reduction(g[i],q)
		f.reverse()
		g.reverse()
		return g

#Subtract g from f in Rq
def poly_sub(f1,g1,q):
	f = f1[:]
	g = g1[:]
	f.reverse()
	g.reverse()
	if len(f)>len(g):
		for i in range(0,len(g)):
			f[i] -= g[i]
			f[i] = modulo_reduction(f[i],q)
		f.reverse()
		g.reverse()
		return f
	else:
		for i in range(0,len(f)):
			g[i] = modulo_reduction((f[i]-g[i]),q)
		for i in range(len(f),len(g)):
			g[i] = modulo_reduction((-g[i]),q)
		f.reverse()
		g.reverse()
		g = remove_lead_zero(g)
		return g

# remove leading zeroes in polynomial
def remove_lead_zero(f):
	while f[0]==0:
		if len(f)==1: return f
		f.remove(0)
	return f

# maps f to a range ([-q/2],[q/2])
def modulo_reduction(f,q):
	f = f%q
	if (f>int(q/2)) : f -=q
	return f

#kind of discrete gaussian distribution, use Modified Bessel function with exp(-t)
def gaussian(n,B_bound):
	j = -1
	while j == -1:
		f = []
		j=0
		for i in range(0,n):
			f.append(modulo_reduction(int((math.exp(-random.randint(0,n))*kn([random.randint(0,10)],1)[0])*10000000),B_bound))
			if f[i]==0 : j+=1
		if j>int(len(f)/2) : j = -1
	f = remove_lead_zero(f)
	return f

#factorisation of a number taken from http://stackoverflow.com/questions/6800193/what-is-the-most-efficient-way-of-finding-all-the-factors-of-a-number-in-python
def factors(n):    
	return set(reduce(list.__add__,([i, n//i] for i in range(1, int(n**0.5) + 1) if n % i == 0)))

#Mobius function taken from https://wiki.math.ntnu.no/_media/ma3001/2014v/analytisktallteori/arithmeticalfunctions_tjerand.pdf
def isValid (n):
	if (n <=0):
		return 0
	elif (n-int (n) != 0):
		return 0
	return 1

def isDivisibleBySquare (n):
	i = 2
	while (i**2 <=n):
		if (n%(i **2)==0):
			return 1
		i += 1
	return 0

def numberOfFactors (n):
	counter = 0
	while (n >1):
		i = 2
		while (i <=n):
			if (n%i ==0):
				n /= i
				counter += 1
				break
			i += 1
	return counter

#The Mobius function .
def mu(n):
	if not isValid(n):
		return 0
	elif (n ==1):
		return 1
	elif (isDivisibleBySquare(n)):
		return 0
	else :
		return ( -1)**numberOfFactors(n)

def mu_poly(d):
	tmp = []
	tmp.append(1)
	for i in range(1,d):
		tmp.append(0)
	tmp.append(-1)
	return tmp
# return cyclotomic polinomial in Rq
def cyclotomicPoly(n,q):
	f = factors(n)
	up = [1]
	down = [1]
	for d in f:
		p_pow = mu(n/d)
		tmp = mu_poly(d)
		if p_pow>0 :
			up = poly_mult(up,tmp,q)
		if p_pow<0 :
			down = poly_mult(down,tmp,q)
	return poly_div(up,down,q)

def check_quo(f,q):
	for i in range(0,len(f)):
		if abs(f[i]) > int(q/4) : return False
	return True

# generate F polynomial
def randomKeyPoly(maxDegree,B_bound):
	output = [] 
	output.append(0)
	for i in range (1,maxDegree):
		random.seed()
		output.append(((int((10**16)*random.random())) % 3)-1)
		#if random.choice('01')=='0' : output[i]=-output[i]
	return output

def take_f_from_cyclo_distribution(n,q):
	f_from_cyclo = random.sample(range(0, n-1), 2)
	f = [1]
	for i in f_from_cyclo:
		f = poly_mult(f,cyclotomicPoly(i,q)[0],q) 
	return poly_add([1],poly_mult([2],f,q),q)

def take_f_from_gaussian_distribution(n,B_bound,q):
	return poly_add([1],poly_mult([2],gaussian(n,B_bound),q),q)
