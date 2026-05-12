from __future__ import annotations

from dataclasses import dataclass
from functools import cached_property, lru_cache
from fractions import Fraction
import math

from nt.rings import ImaginaryQuadraticField, QuaternionAlgebra
from nt.common import (
    legendre,
    factorize,
    elements_of_exact_order,
    valuation as vl,
    divisors,
    fmt_factored,
    fmt_invariants,
    sgn,
    equiv,
    coprime_part,
)

from utils.logging import Logger

@dataclass
class ConductorRecord:
    """Per-conductor contribution within one isogeny class at level N."""
    f: int
    n_curves: int            # h(f), unweighted curve count
    n_pts_exact_order: int   # |{x in E(F_q)[N] : ord(x) = N}|
    value: Fraction          # n_pts_exact_order * hw(f), the weighted contribution


@dataclass
class TorsionRecord:
    """Aggregate torsion data for one isogeny class at level N."""
    n_curves: int
    n_points: int
    value: Fraction
    conductor_levels: list[ConductorRecord]  # per-conductor breakdown (for debugging)


class IsogenyClass:
    def __init__(self, p: int, n:int, t: int, N_t: int = None) -> None:
        self.q = p**n
        self.p = p
        self.n = n

        self.t = t
        self.a_pi = (self.t-self.f_pi) // 2 if self.field.D % 4 == 1 else (self.t // 2)

        self.__N_t = N_t
        self.__torsion_records: dict[int, TorsionRecord] = {}   # keyed by N
        self.__group_cache: dict[int, tuple[int, int]] = {}     # f -> full (e1, e2)

    @cached_property
    def D_pi(self) -> int:
        return self.t**2 - 4 * self.q

    @cached_property
    def f_pi(self) -> int:
        return math.isqrt(abs(self.D_pi // self.field.D)) if self.field.D != 0 else 1

    @cached_property
    def ordinary(self) -> bool:
        return self.t % self.p != 0

    @cached_property
    def n_pts(self) -> int:
        return self.q + 1 - self.t

    @cached_property
    def n_pts_total(self) -> int:
        return self.n_pts * self.size

    @cached_property
    def field(self) -> ImaginaryQuadraticField:
        return (
            ImaginaryQuadraticField(self.D_pi) if self.D_pi < 0 else QuaternionAlgebra(self.p)
        )

    @cached_property
    def is_quaternion(self) -> bool:
        return self.D_pi == 0

    @cached_property
    def size(self) -> int:
        return (
            self.__N_t
            if self.__N_t is not None
            else (ImaginaryQuadraticField.H(self.D_pi) if self.D_pi < 0 else QuaternionAlgebra.H(self.p))
        )

    @lru_cache(maxsize=4096)
    def eval_hk(self, k: int) -> int:
        return int(sum(
            math.comb(k - j, j) * (-self.q) ** j * self.t ** (k - 2 * j)
            for j in range(k // 2 + 1)
        ))

    def eval_torsion(self, N: int, flatten: bool = True) -> Fraction:
        return self._torsion_record(N, flatten).value

    def get_torsion(self, N: int, flatten: bool = True) -> TorsionRecord:
        return self._torsion_record(N, flatten)

    def _torsion_record(self, N: int, flatten: bool = True) -> TorsionRecord:
        if N not in self.__torsion_records:
            self.__torsion_records[N] = self._compute_torsion_record(N, flatten)
        return self.__torsion_records[N]

    def _l_sylow(self, l: int, f: int) -> tuple[int, int]:
        hf = max(0, vl(self.f_pi, l) - vl(f, l))
        r1 = min(hf, vl(self.n_pts, l) // 2)  # hf was listed twice before
        return (r1, vl(self.n_pts, l) - r1)

    def _full_group_structure(self, f: int) -> tuple[int, int]:
        """Full invariants (e1, e2) of E(F_q) for curves with conductor f."""
        if f not in self.__group_cache:
            if self.is_quaternion:
                e1 = self.p**(self.n // 2) - sgn(self.t)
                e2 = e1
            else:
                e1 = math.gcd(self.a_pi-1, self.f_pi // f)
                e2 = self.n_pts // e1
                '''e1, e2 = 1, 1
                for p, _ in factorize(self.n_pts):
                    r1, r2 = self._l_sylow(p, f)
                    e1 *= p ** r1
                    e2 *= p ** r2'''
            self.__group_cache[f] = (e1, e2)
        return self.__group_cache[f]

    def _torsion_subgroup(self, f: int, N: int, compute_full:bool = False) -> tuple[int, int]:
        """N-torsion subgroup invariants for curves with conductor f.

        Uses the cached full group structure if available (gcd path),
        otherwise computes restricted to the prime factors of N only.
        """
        if f in self.__group_cache:
            e1, e2 = self.__group_cache[f]
        else:
            e1, e2 = self._full_group_structure(f)  
            # elif self.is_quaternion or compute_full:
            #    e1, e2 = self._full_group_structure(f) # we can not use conductors for the quaternion
            # _e1 = math.gcd(math.gcd(self.a_pi-1, self.f_pi // f), N )
            # _e2 = math.gcd(self.n_pts, N // _e1)

            '''e1, e2 = 1, 1
            for p, _ in factorize(N):
                r1, r2 = self._l_sylow(p, f)
                e1 *= p ** r1
                e2 *= p ** r2'''
            # inv1 = fmt_invariants((math.gcd(N, e1), math.gcd(N, e2)))
            # inv_ = fmt_invariants((_e1, _e2))
            # grp_inv = self._full_group_structure(f)
            # fmt_invariants(grp_inv)

        return (math.gcd(N, e1), math.gcd(N, e2))

    def _compute_torsion_at(self, f: int, N: int, compute_full: bool = False) -> ConductorRecord:
        torsion_inv = self._torsion_subgroup(f, N, compute_full)

        grp_inv = self._full_group_structure(f)
        p_inert = (
            2
            if (
                not self.ordinary
                and not self.is_quaternion
                and legendre(self.field.D * f**2, self.p) < 0
            )
            else 1
        )
        n_pts = elements_of_exact_order(N, torsion_inv[0], torsion_inv[1])
        n_curves = self.field.h(f) * p_inert
        n_curves_weighted = self.field.hw(f) * p_inert

        return ConductorRecord(
            f=f,
            n_curves=n_curves,
            n_pts_exact_order=n_pts,
            value=n_pts * n_curves_weighted,
        )

    def _compute_torsion_record(self, N: int, flatten: bool = True) -> TorsionRecord:
        if self.n_pts % N != 0:
            return TorsionRecord(n_curves=0, n_points=0, value=Fraction(0), conductor_levels=[])

        # Strip p-part unconditionally — conductors divisible by p are always excluded
        # TODO: need to double check this applies to not only SS
        f_pi_reduced = self.f_pi // self.p ** vl(self.f_pi, self.p)
        H_coprime = 1

        # Flatten: collapse the N-coprime part of f_pi into a scalar H_coprime.
        # Only valid for generic imaginary quadratic fields (D_K < -4) — Eisenstein
        # and Gaussian fields have non-standard aut groups at f=1 that break
        # multiplicativity of h(O_f) in the coprime tower.
        flatten = True
        # p = 17 test is good

        clr = Logger.HEADLINE
        if self.f_pi % N == 0:
            clr = Logger.FAIL if self.n_pts % (N * N) != 0 else Logger.NOTICE
        '''Logger.cprint(
            f"Isogeny class with t={self.t}, f_pi={fmt_factored(self.f_pi)}, n_pts={fmt_factored(self.n_pts)}, N={N}, q equiv N={equiv(self.q,N)}",
            clr
        )'''

        can_flatten = flatten and not self.is_quaternion
        padding = 1

        if can_flatten:
            coprime = coprime_part(f_pi_reduced, N)
            f_pi_reduced //= coprime
            f_list = divisors(f_pi_reduced)
            '''print(
                f"Flattening conductor from {fmt_factored(self.f_pi)} to {fmt_factored(f_pi_reduced)}, coprime={fmt_factored(coprime)}, f_list={f_list}"
            )'''
            if coprime > 1:
                H_coprime = self.field.Hf_inv(coprime)
        else:
            f_list = divisors(f_pi_reduced)
        conductor_levels = [
            self._compute_torsion_at(f, N) for f in f_list
        ]

        return TorsionRecord(
            n_curves= H_coprime * sum(c.n_curves for c in conductor_levels),
            n_points= H_coprime * sum(c.n_curves * c.n_pts_exact_order for c in conductor_levels),
            value= H_coprime * sum(c.value for c in conductor_levels),
            conductor_levels=conductor_levels,
        )


class CurvesRecordFq:
    """All isogeny classes of elliptic curves over F_q."""

    def __init__(self, p: int, n: int) -> None:
        self.q = p**n
        self.p = p
        self.n = n
        self.j0_SS = legendre(-3, p) == -1  # note: true for p = 2, not true for p = 3
        self.j1728_SS = (
            legendre(-4, p) == -1
        )  # note: true for p = 2, not true for p = 3

        self.__isogeny_classes: dict[int, IsogenyClass] = {
            t: IsogenyClass(self.p, self.n, t, N_t)
            for t, N_t in self.__ordinary_ts + self.__supersingular_ts
        }

    @cached_property
    def HB(self) -> int:
        return math.isqrt(4 * self.q)

    @cached_property
    def __ordinary_ts(self) -> list[(int, int)]:
        """
        Ordinary traces t satisfy |t| < 2*sqrt(q) and p does not divide t.
        N(t) = H(t^2 - 4q) counts the number of curves in these isogeny classes, no inert factor scaling for ordinary curves.
        see: Schoof 1987, Thm 4.6 for reference

        Return a list of tuples (t, None) for ordinary traces, where N(t) will fallback to H(t^2 - 4q). NOTE: this is NOT weighted by 1/|Aut(E)|

        """

        return [
            (s * t, None)
            for t in range(1, self.HB + 1)
            if t % self.p != 0
            for s in (1, -1)
        ]

    @cached_property
    def __supersingular_ts(self) -> list[(int, int)]:
        """
        Supersingular traces t satisfy specific conditions depending on the characteristic and degree of the field.
        see: Schoof 1987, Thm 4.6 for reference on each case

        Return a list of tuples (t, N(t)) for supersingular traces, where N(t) is the number of curves in the isogeny class with trace t, including any inert factor scaling. NOTE: this is NOT weighted by 1/|Aut(E)|
        """

        ts: list[(int, int)] = []
        if self.n % 2 == 1:
            """
            For char(Fq) > 3, the only permitted SS trace is 0.
            This IsogenyClass will lie in field of D_K = -4p, and all curves will lie in max order, hence: N(t) = H(-4p)
            """
            ts.append((0, ImaginaryQuadraticField.H(-4 * self.p)))
            """
            For char(Fq) = 2, 3, we get t = pm sqrt(2q) or t = pm sqrt(3q) respectively, producing exactly one curve per such trace, N(t) = 1.
            """
            if self.p == 2 or self.p == 3:
                t = self.p ** ((self.n + 1) // 2)
                ts += [(s * t, 1) for s in (1, -1)]
        else:

            # N(t) = 1 - (-3/p), will yield D_K = -3 and only max order allowed, only curve here is j0 with inert factor 2
            if self.j0_SS:
                sqrt_q = self.p ** (self.n // 2)
                ts += [(s * sqrt_q, 2) for s in (1, -1)]

            # N(t) = 1 - (-4/p), will yield D_K = -4 and only max order allowed, only curve here is j1728 with inert factor 2
            if self.j1728_SS:
                ts.append((0, 2))

            # Quaternions, we compute the full class number as given in Schoof 1987, Thm 4.6
            ts += [(s * self.HB, QuaternionAlgebra.H(self.p)) for s in (1, -1)]
        return ts

    @cached_property
    def num_curves_total(self) -> Fraction:
        return Fraction(0) + sum(c.N_t for c in self.isogeny_classes.values())

    @cached_property
    def isogeny_classes(self) -> Fraction:
        return self.__isogeny_classes.items()

    def check(self) -> bool:
        return self.num_curves_total == self.q
