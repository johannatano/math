from __future__ import annotations
from abc import ABC, abstractmethod
from functools import cached_property
from typing import Sequence


class AbelianGroup(ABC):
    """Abstract finite abelian group G ≅ Z/n1 x Z/n2 (with n1 | n2)."""

    @property
    @abstractmethod
    def invariants(self) -> tuple[int, int]:
        """Return (n1, n2) with n1 | n2 such that G ≅ Z/n1 x Z/n2."""

    @abstractmethod
    def order(self) -> int:
        """Return |G|."""

    @abstractmethod
    def points_of_exact_order(self, n: int) -> int:
        """Return the number of elements of exact order n in G."""

    @abstractmethod
    def p_sylow(self, p: int) -> AbelianGroup:
        """Return the p-Sylow subgroup of G."""


class AbelianGroup_ZnxZm(AbelianGroup):
    """Concrete finite abelian group given explicitly by its invariants (n1, n2)."""

    def __init__(self, n1: int, n2: int) -> None:
        assert n1 == 1 or n2 % n1 == 0, "require n1 | n2"
        self._n1 = n1
        self._n2 = n2

    @property
    def invariants(self) -> tuple[int, int]:
        return (self._n1, self._n2)

    def order(self) -> int:
        return self._n1 * self._n2

    def points_of_exact_order(self, n: int) -> int:
        from math import gcd
        # |{x in Z/n1 x Z/n2 : ord(x) = n}| = gcd(n,n1)*gcd(n,n2) - gcd(n/p,n1)*gcd(n/p,n2) summed...
        # use inclusion-exclusion over prime divisors of n
        from sympy import factorint
        def phi_group(k: int) -> int:
            """Elements whose order divides k."""
            return gcd(k, self._n1) * gcd(k, self._n2)
        result = 0
        # Möbius inversion: #{ord = n} = sum_{d|n} μ(n/d) * #{ord | d}
        from sympy import mobius, divisors
        for d in divisors(n):
            result += mobius(n // d) * phi_group(d)
        return result

    def p_sylow(self, p: int) -> AbelianGroup_ZnxZm:
        def p_part(k: int) -> int:
            while k % p == 0:
                k //= p
            return self._n1 * self._n2 // k  # p-power part
        from math import gcd
        # p-Sylow of Z/n1 x Z/n2 is Z/p^a1 x Z/p^a2
        def vp(k: int) -> int:
            v = 0
            while k % p == 0:
                k //= p
                v += 1
            return v
        a1, a2 = vp(self._n1), vp(self._n2)
        return AbelianGroup_ZnxZm(p**a1, p**a2)
