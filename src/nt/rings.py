from __future__ import annotations
from dataclasses import dataclass

from functools import cached_property, lru_cache
from fractions import Fraction
import math
from sage.libs.pari import pari

from pathlib import Path
from nt.common import legendre, factorize, divisors

_DB_PATH = Path(__file__).parent / "store/classnb.db"

# Set True to save pari-computed values back to the db on the fly.
_DB_SAVE = False

# Lazily opened read connection — None until first db lookup.
_db_conn = None
max_D = -1

def _db_lookup(D: int):
    """Return h(D) from db, or None if not found / db absent."""
    global _db_conn
    if not _DB_PATH.exists():
        return None
    if _db_conn is None:
        import sqlite3
        _db_conn = sqlite3.connect(f"file:{_DB_PATH}?mode=ro", uri=True)
    try:
        row = _db_conn.execute("SELECT h FROM classnb WHERE D=?", (D,)).fetchone()
        return row[0] if row else None
    except Exception:
        return None


@lru_cache(maxsize=None)
def _H(D: int) -> int:
    global max_D
    if abs(D) > max_D:
        max_D = abs(D)
        h = _db_lookup(D)
    h = _db_lookup(D)
    if h is not None:
        return h
    h = int(pari(D).qfbclassno())
    if _DB_SAVE:
        import sqlite3
        conn = sqlite3.connect(_DB_PATH)
        conn.execute("INSERT OR IGNORE INTO classnb VALUES (?,?)", (D, h))
        conn.commit()
        conn.close()
    return h


class Order:
    def __init__(self, K: NumberField, f: int):
        self.K = K
        self.f = f

    @property
    def is_maximal(self) -> bool:
        return self.f == 1

    @cached_property
    def D(self) -> int:
        return self.f**2 * self.K.D

    @cached_property
    def unit_index(self) -> int:
        """[O_K^× : O_f^×] — equals 1 for f>1 except field has extra roots of unity."""
        # if self.is_maximal:
        #    return 1
        return self.K._maximal_unit_index

    @cached_property
    def h(self) -> int:
        raise NotImplementedError

    @cached_property
    def hw(self) -> int:
        raise NotImplementedError

    def order_of(self, f: int) -> Order:
        """Return the sub-order of conductor f (must divide self.f)."""
        return self.K.order(f)

class QuadraticOrder(Order):
    @cached_property
    def h_prime_prod(self) -> Fraction:
        ##if self.is_maximal:
        ##    return Fraction(1)
        result = Fraction(self.f)
        for p, _ in factorize(self.f):
            result *= Fraction(p - legendre(self.K.D, p), p)
        return result

    @cached_property
    def h(self) -> int:
        if self.is_maximal:
            return _H(self.K.D)
        '''print(
            f"Computing h={Fraction(self.K.O_K.h, self.unit_index) * self.h_prime_prod}"
        )'''
        return Fraction(self.K.O_K.h, 1) * self.h_prime_prod
        ##return Fraction(self.K.O_K.h, 1) * self.h_prime_prod
        '''for p, _ in factorize(self.f):
            result *= Fraction(p - legendre(self.K.D, p), p)
        assert (
            result.denominator == 1
        ), f"Non-integer class number: D_K={self.K.D}, f={self.f}, result={result}"
        return int(result)'''

    @cached_property
    def hw(self) -> int:
        """Returns weight by 1/|Aut| (Schoof 1987, Thm 4.6)."""
        # self.K._maximal_unit_index * 2 if self.is_maximal else
        return Fraction(
            self.h,
            self.K._maximal_unit_index * 2,  # this re-normalizes the class number into weighted by aut
        )


class QuaternionOrder(Order):
    @cached_property
    def h(self) -> int:
        """Curve count in Q_{∞,p} (Schoof 1987, Thm 4.6)."""
        if self.f != 1:
            raise NotImplementedError("Non-maximal quaternion orders not implemented")
        return Fraction(
            self.K.p + 6 - 4 * legendre(-3, self.K.p) - 3 * legendre(-4, self.K.p), 12
        )

    @cached_property
    def hw(self) -> int:
        """Curve count in Q_{∞,p} weighted with 1/|Aut| (Schoof 1987, Thm 4.6)."""
        if self.f != 1:
            raise NotImplementedError("Non-maximal quaternion orders not implemented")
        return Fraction(self.K.p - 1, 24)


class NumberField:
    @staticmethod
    def _D0(D: int) -> tuple[int, int]:
        """Return (D_K, f) where D_K is the fundamental discriminant and D = D_K * f²."""
        if D == 0:
            return 0, 1
        sign = -1 if D < 0 else 1
        from nt.common import factorize
        sf = 1
        for p, e in factorize(abs(D)):
            if e % 2 == 1:
                sf *= p
        sf *= sign
        D_K = sf if sf % 4 == 1 else 4 * sf
        f = math.isqrt(abs(D) // abs(D_K))
        return D_K, f

    def __init__(self, D: int) -> None:
        self._disc, self._conductor = NumberField._D0(D)
        self._order_cache: dict[int, Order] = {}
        self._maximal_unit_index = 1

    @cached_property
    def D(self) -> int:
        return self._disc

    @cached_property
    def O_K(self) -> Order:
        return self.order(1)

    def h(self, f:int) -> int:
        return self.order(f).h

    def h_tilde(self, f: int) -> int:
        return self.order(f).h

    def hw(self, f: int) -> int:
        return self.order(f).hw

    def Hf(self, f: int) -> int:
        raise NotImplementedError

    def _make_order(self, f: int) -> Order:
        raise NotImplementedError

    def order(self, f: int) -> Order:
        if f not in self._order_cache:
            self._order_cache[f] = self._make_order(f)
        return self._order_cache[f]

    @cached_property
    def u_index(self) -> int:
        return self._maximal_unit_index


# --- Imaginary quadratic fields and quaternion algebras ---


class ImaginaryQuadraticField(NumberField):

    @staticmethod
    @lru_cache(maxsize=4096)
    def H2(D: int) -> int:
        """Kronecker class number H(Δ) — unweighted count of SL2(Z)-orbits (Schoof 1987)."""
        if D >= 0:
            return 0
        total = 0
        abs_D = abs(D)
        for d in range(1, math.isqrt(abs_D) + 1):
            if abs_D % (d * d) == 0:
                delta = D // (d * d)
                if delta % 4 == 0 or delta % 4 == 1:
                    total += _H(delta)
        return total

    @staticmethod
    @lru_cache(maxsize=4096)
    def H(D: int) -> int:
        # to what fundamental discriminant do we belong, and at what height do we probe
        DO, f0 = NumberField._D0(D)
        # just return the cached class number of the order at that height
        return ImaginaryQuadraticField(DO).Hf(f0)

    def __init__(self, D: int) -> None:
        super().__init__(D)
        self._maximal_unit_index = (
            3 if self._disc == -3 else (2 if self._disc == -4 else 1)
        )

    def Hf(self, f: int) -> int:

        for e in divisors(f):
            print(
                f"  divisor e={e}, h(e)={self.order(e).h}, unit_index={self.order(e).unit_index}, adding={self.order(e).h}"
        )
        return sum(self.order(e).h for e in divisors(f))

    def Hf_inv(self, f: int) -> int:
        """H(D_K · f²) using cached order class numbers — avoids pari calls.
        H'(D_K · f²)  = (Σ_{e | f} h(e)*h(e).unit_idx) // h(1)
        """
        return sum(self.order(e).h for e in divisors(f)) // self.O_K.h

    def Hf_tilde(self, f: int) -> int:
        """H(D_K · f²) using cached order class numbers — avoids pari calls.
        H'(D_K · f²)  = (Σ_{e | f} h(e)*h(e).unit_idx) // h(1)
        """
        return sum(self.order(e).h_prime_prod for e in divisors(f))

    def _make_order(self, f: int) -> QuadraticOrder:
        return QuadraticOrder(self, f)

    def h_tilde(self, f: int) -> int:
        return self.order(f).h_prime_prod


class QuaternionAlgebra(NumberField):
    @staticmethod
    @lru_cache(maxsize=4096)
    def H(p:int, normalized: bool = False) -> Fraction:
        """Curve count in Q_{∞,p} (Schoof 1987, Thm 4.6).
        normalized=False: weighted count with 1/|Aut| for j=0 and j=1728.
        normalized=True:  full mass, simplifies to (p-1)/12.
        """
        if normalized:
            return Fraction(p - 1, 12)
        chi3 = legendre(-3, p)
        chi4 = legendre(-4, p)
        return Fraction(p + 6 - 4 * chi3 - 3 * chi4, 12)

    def __init__(self, p: int) -> None:
        super().__init__(0)
        self.p = p

    def Hf(self, f: int) -> QuaternionOrder:
        return QuaternionOrder.H(self.p)
    

    def _make_order(self, f: int) -> QuaternionOrder:
        return QuaternionOrder(self, f)
