q = 0x1000000000000000000000000000000014def9dea2f79cd65812631a5cf5d3ed
h = 10

# c = 2^550
c = 3685510180489786476798393145496356338786055879312930105836138965083617346086082863365358130056307390177215209990980317284932211552660930305235775636164742230126366470
n = c * q^h
while (n+1 not in Primes()):
    c = c + 1
    n = c * q^h

p = n+1
Zp = GF(p)
print(p)

while true:
    g = Zp.random_element()^c
    flag = true
    for i in range(1, h):
        print(i)
        if pow(g, q^i) == 1:
            flag = false
    if (flag and pow(g, q^h)) == 1:
        break

def elgamal_enc(m, pk):
    r = randint(1, q^h)
    c1 = g^r
    c2 = g^m * pk^r
    return (c1, c2)

def elgamal_dec(c1, c2, sk):
    m = c2 * c1^-sk
    return m

sk = randint(1, q^h)
pk = g^sk

(c1, c2) = elgamal_enc(2, pk)
(d1, d2) = elgamal_enc(3, pk)

print(elgamal_dec(c1, c2, sk)) == g^2
print(elgamal_dec(d1, d2, sk)) == g^3

c1 *= d1
c2 *= d2

print(elgamal_dec(c1, c2, sk)) == g^5

print(c1^(q^h))
print(c2^(q^h))

