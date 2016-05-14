import math

#Funciones de probabilidad
#1
def hiper(nn, m, n, r, flag):# Si flag es 1 se calcula probabilidad puntual, si es flag 0 acumulada
    return sum(float(comb(m, x) * comb(nn-m, n-x)) / comb(nn, n) for x in range(r * flag,r + 1))
#2
def binom(n, p, r, flag):# Si flag es 1 se calcula prbabilidad puntual, si es flag 0 acumulada
    return sum(comb(n, x) * (p ** x) * ((1-p)**(n-x)) for x in range(r * flag,r + 1))
#3
def binom_negativa(r, k, p, flag):
    return sum(comb(x - 1, k - 1) * (p ** k) * ((1 - p) ** (x - k)) for x in range(r * flag,r + 1))
#4
def poisson(r, lammbda, flag):
    return sum(float((lammbda ** x) * math.exp(-lammbda)) / fact(x) for x in range(r * flag,r + 1))
#5
def exponencial(r, lammbda, flag):
    return sum(1 - (math.e ** (-lammbda * t)) for t in range(r * flag,r + 1))
#6
def z_normal(d, m, x):
    return (x-m)/d
#7
def norm_hastings(z, flag):
    if (flag == 1):# Si flag es 1 calcula la puntual
        return F(z)
    else:#Si flag es 0 calcula la acumulada
        p0 = 0.2316419
        c1 = 0.31938153
        c2 = -0.356563781
        c3 = 1.781477937
        c4 = -1.821255979
        c5 = 1.330274429

        if (z >= 0):
            t = 1 / (1 + p0 * z)
            return 1 - F(z) * (c1*t + c2*t**2 + c3*t**3 + c4*t**4 + c5*t**5)
        else: #Si z es menor que 0 se multiplica por -1 para calcular t
            t = 1 / (1 + p0 * -z)
            return F(z) * (c1*t + c2*t**2 + c3*t**3 + c4*t**4 + c5*t**5)                                     
#8
def binom_normal(p,n,x):
    if (p * n >= 5) and ((1-p) * n >= 5):
        d = math.sqrt(n * p * (1-p))
        z = (x - n*p)/d
        return norm_hastings(z, 0)
    print "No se puede aproximar usando normal"
#9
def triangular(x, a, b, c, flag):
    if flag == 1:
        if x < a:
            return 0
        if x >= a and x <= c:
            return float(2 * (x-a))/((b-a)*(c-a))
        if x > c and x <= b:
            return float(2 * (b - x))/((b-a)*(b-c))
        return 0
    else:
        if x < a:
            return 0
        if x >= a and x <= c:
            return float((x-a)**2)/((b-a)*(c-a))
        if x > c and x < b:
            return 1 - (float((b-x)**2)/((b-a)*(b-c)))
        return 1
#10                           
def uniforme(a, b, x, flag):
    if flag == 1:
        if x <= b and x >= a:
            return 1/float(b-a)
        return 0
    else:
        if x < a:
            return 0
        if x <= b and x >= a:
            return (x-a)/(b-a)
        return 1

#Funciones inversas
#11
def inv_exponencial(u, lammbda):
    return math.log(1 - u)/-lammbda
#12
def inv_uniforme(a, b, u):
    return (b-a)*u + a
#13
def inv_triangular(u, a, b, c):
    if u >= 0 and u <= (c-a)/(b-a):
        return math.sqrt(u*(b-a)*(c-a)) + a
    if u > (c-a)/(b-a) and u <= 1:
        return b - math.sqrt((1-u)*(b-a)*(b-c))
#14
def inv_normal(p):#funcion probit
    return math.sqrt(2) * inv_erf(2*p - 1)

#Funciones auxiliares
#15
def fact(num):
    resp = 1
    for x in range(2, num + 1):
    	resp = resp * x
    return resp
#16
def comb(n, r):
    if r > n or r < 0 or n < 0:
        return 0
    else:
        return fact(n) / (fact(n - r) * fact(r))
#17
def F(z):
    return math.exp(-z**2 / 2) / math.sqrt(2 * math.pi)
#18
def inv_erf(v):
    x = (math.sqrt(math.pi)*v)/2
    return (x + (1.0/3)*(x**3) + (7.0/30)*(x**5) + (127.0/630)*(x**7) + (4369.0/22680)*(x**9) + (34807.0/178200)*(x**11) + (20036983.0/97297200)*(x**13))

print "1.- hiper(80, 5, 6, 0, 1) = " + str(hiper(80, 5, 6, 0, 1))
print "2.- binom(6, 0.0625, 0, 1) = " + str(binom(6, 0.0625, 0, 1))
print "3.- binom_negativa(6, 3, 0.05, 1) = " + str(binom_negativa(6, 3, 0.05, 1))
print "4.- poisson(3, 1.1934, 1) = " + str(poisson(3, 1.1934, 1))
print "5.- exponencial(1, 2.1, 1) = " + str(exponencial(1, 2.1, 1))
print "6.- z_normal(1.3, 45, 42) = " + str(z_normal(1.3, 45, 42))
print "7.- norm_hastings(-2.307, 0) = " + str(norm_hastings(-2.307, 0))
print "8.- binom_normal(0.5,600,249) = " + str(binom_normal(0.5,600,249))
print "9.- triangular(3, 2, 7, 5, 1) = " + str(triangular(3, 2, 7, 5, 1))
print "10.- uniforme(2, 7, 3, 1) = " + str(uniforme(2, 7, 3, 1))
print "11.- inv_exponencial(0.056, 3) = " + str(inv_exponencial(0.056, 3))
print "12.- inv_uniforme(5, 7,0.355) = " + str(inv_uniforme(5, 7,0.355))
print "13.- inv_triangular(0.5, 2, 7, 5) = " + str(inv_triangular(0.5, 2, 7, 5))
print "14.- inv_normal(0.0107) = " + str(inv_normal(0.0107))
print "15.- fact(4) = " + str(fact(4))
print "16.- comb(5, 4) = " + str(comb(5, 4))
print "17.- F(-2.307) = " + str(F(-2.307))
print "18.- inv_erf(0.5) = " + str(inv_erf(0.5))
