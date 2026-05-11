from __future__ import annotations
from functools import lru_cache

def euler_phi(n):
    result = n
    for p, _ in factorize(n):
        result -= result // p
    return result


def legendre(a, p):
    if p == 2:
        return 0 if a % 2 == 0 else (1 if a % 8 in (1, 7) else -1)
    a = a % p
    if a == 0:
        return 0
    r = pow(a, (p - 1) // 2, p)
    return -1 if r == p - 1 else 1


@lru_cache(maxsize=None)
def valuation(n, l):
    if n == 0:
        return float("inf")
    v = 0
    while n % l == 0:
        n //= l
        v += 1
    return v


def fmt_invariants(inv: tuple) -> str:
    a, b = fmt_factored(inv[0]), fmt_factored(inv[1])
    return f"({a} \033[33mx\033[0m {b})"


def fmt_factored(n):
    if n == 0:
        return "0"
    if n > 10**20:
        return str("")
    factors = factorize(abs(n))
    if not factors:
        return str(n)
    s = "·".join(f"{p}^{e}" if e > 1 else str(p) for p, e in factors)
    return f"-{s}" if n < 0 else s


@lru_cache(maxsize=None)
def factorize(n):
    # return list(factorint(n).items())  # → [(2, 2), (3, 1)]
    """Return list of (prime, exponent) pairs for n > 1."""
    factors = []
    d = 2
    while d * d <= n:
        if n % d == 0:
            e = 0
            while n % d == 0:
                n //= d
                e += 1
            factors.append((d, e))
        d += 1
    if n > 1:
        factors.append((n, 1))
    return factors


def divisors(n):
    divs = []
    i = 1
    while i * i <= n:
        if n % i == 0:
            divs.append(i)
            if i != n // i:
                divs.append(n // i)
        i += 1
    return sorted(divs)


def elements_of_exact_order_exp(ell: int, a: int, e1: int, e2: int) -> int:
    """Number of elements of exact order ell^a in Z/ell^e1 x Z/ell^e2."""
    s1 = min(a, e1)
    s2 = min(a, e2)
    return ell ** (s1 + s2) - ell ** (min(a - 1, s1) + min(a - 1, s2))


def elements_of_exact_order(N: int, n1: int, n2: int) -> int:
    """Number of elements of exact order N in Z/n1 x Z/n2.

    Factors over primes: at each ell, count elements whose ell-Sylow part
    has order exactly ell^a. These conditions are independent across primes.
    """
    result = 1
    for ell, a in factorize(N):
        result *= elements_of_exact_order_exp(ell, a, valuation(n1, ell), valuation(n2, ell))
    return result


def sgn(x):
    return 1 if x > 0 else -1 if x < 0 else 0
