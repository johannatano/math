# external libraries
from dataclasses import dataclass
from functools import cached_property
from fractions import Fraction
import math

from sage.libs.pari import pari

# custom libraries
from nt.functions import factorize,legendre

# HELPER METHODS RELATED TO QUADRATIC FIELDS


class QuaternionOrder:
    @staticmethod
    def class_number(p, rescale_weights=False):
        """print(f"p={p}, chi3={chi3}, chi4={chi4}")
        H_p = (p + 6 - 4 * chi3 - 3 * chi4) / 12
        print((p + 6 - 4 * chi3 - 3 * chi4))
        return Fraction(int(p + 6 - 4 * chi3 - 3 * chi4), 12)
        n_j0 = (1 - chi3) // 2
        n_j1728 = (1 - chi4) // 2
        n_generic = H_p - n_j0 - n_j1728"""
        # just return N(t)
        if not rescale_weights:
            chi3 = legendre(-3, p)
            chi4 = legendre(-4, p)
            return Fraction(int(p + 6 - 4 * chi3 - 3 * chi4), 12)
        return Fraction(int(p - 1), 12)


'''
def quadratic_order_class_number(disc):
    r"""
    Return the class number of the quadratic order of given discriminant.

    EXAMPLES::

        sage: from sage.rings.number_field.order import quadratic_order_class_number
        sage: quadratic_order_class_number(-419)
        9
        sage: quadratic_order_class_number(60)
        2

    ALGORITHM: Either :pari:`qfbclassno` or :pari:`quadclassunit`,
    depending on the size of the discriminant.
    """
    # cutoffs from PARI documentation
    if disc < -10**25 or disc > 10**10:
        h = pari.quadclassunit(disc)[0]
    else:
        h = pari.qfbclassno(disc)
    return ZZ(h)'''

@dataclass
class Order:
    f: int
    D: int
    h_ord: int  # class nr for this order
    inv: Tuple[int, int]  # all curves in this order have same invariants


class QuadraticField:
    def __init__(self, D: int) -> None:
        self.D_K = self._fundamental_discr(D)
        self.h_OK = -1 if self.D_K == 0 else int(pari(self.D_K).qfbclassno())
        # Unit index [O_K^x : O_f^x]
        if self.D_K == -3:
            self.w = 3
        elif self.D_K == -4:
            self.w = 2
        else:
            self.w = 1
    def _fundamental_discr(self, D: int) -> int:
        if D == 0:
            return 0
        sign = -1 if D < 0 else 1
        abs_D = abs(D)
        # Remove all square factors to find squarefree part
        sf = 1
        n = abs_D
        d = 2
        while d * d <= n:
            exp = 0
            while n % d == 0:
                n //= d
                exp += 1
            if exp % 2 == 1:
                sf *= d
            d += 1
        if n > 1:
            sf *= n
        sf *= sign  # restore sign
        # Fundamental discriminant: sf if sf ≡ 1 mod 4, else 4·sf
        delta = sf if sf % 4 == 1 else 4 * sf
        # f = math.isqrt(abs(D // delta))
        return int(delta)  # , f

    def h(self, f: int) -> int:
        """
        Class number of the order of conductor f in the imaginary quadratic field
        with fundamental discriminant D_K.
        Matches Sage's K.order_of_conductor(f).class_number().
        Parameters:
            f    : conductor (positive integer)
        """
        if f == 1:
            return self.h_OK
        # Product over distinct primes dividing f
        
        result = f * self.h_OK // self.w
        
        for p, _ in factorize(f):
            result = result * (p - legendre(self.D_K, p) // p)
        return int(result)

    def j_invariants(self, f: int, F) -> List:
        j_invs = []
        """try:
            H = hilbert_class_polynomial(self.D_K*f**2)
            H_fq = H.change_ring(F)
            for j, m in H_fq.roots(multiplicities=True):
                for _ in range(m):
                    j_invs.append(j)
        except Exception as e:
            print(f"Warning: Could not compute HCP for D={self.D_K*f**2}, probably because sage is not loaded: {e}")"""
        return j_invs


@dataclass
class IsogenyClassLevelResult:
    level: int
    ''' note these can be fractions due to 1 / #aut weight '''
    total_num_curves: Fraction = field(default_factory=Fraction) 
    total_num_SS: Fraction = field(default_factory=Fraction)
    contrib_value: Fraction = field(default_factory=Fraction)
    ''' This will always be integer '''
    total_num_points: int

@dataclass
class IsogenyClassOrder:
    f : int
    inv: Tuple[int, int] # all curves in here will have same invariants
    num_curves: int

class IsogenyClass:
    def __init__(self, t: int, q:int) -> None:
        self.t = t
        self.q = q
        self.D_pi = t**2 - 4 * q
        self.N_pts = q + 1 - t # E
        self.N = 0 # N(t) = H(D_pi)
    
    # cache
    def hk(self, k: int) -> int:
        return int(sum(
            math.comb(k - j, j) * (-self.q) ** j * self.t ** (k - 2 * j)
            for j in range(k // 2 + 1)
        ))
