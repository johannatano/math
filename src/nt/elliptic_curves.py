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
)

from utils.logging import Logger

@dataclass
class TorsionRecord:
    n_curves: int
    n_points: int
    value: Fraction # this is the normalized value total number fo points / aut


class IsogenyClass:
    def __init__(self, p: int, n:int, t: int, N_t: int = None) -> None:
        self.q = p**n
        self.p = p
        self.n = n

        self.t = t
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
            else ImaginaryQuadraticField.H(self.D_pi)
        )

    @lru_cache(maxsize=4096)
    def eval_hk(self, k: int) -> int:
        return int(sum(
            math.comb(k - j, j) * (-self.q) ** j * self.t ** (k - 2 * j)
            for j in range(k // 2 + 1)
        ))

    @lru_cache(maxsize=4096)
    def eval_torsion(self, N: int) -> int:
        if N not in self.__torsion_records:
            self.__torsion_records[N] = self._compute_torsion_record(N)
        return self.__torsion_records[N].value

    @lru_cache(maxsize=4096)
    def get_torsion(self, N: int) -> int:
        if N not in self.__torsion_records:
            self.__torsion_records[N] = self._compute_torsion_record(N)
        return self.__torsion_records[N]

    def _l_sylow(self, l:int, f:int) -> dict[int, int]:
        # NOTE: this returns the exponents only
        hf = max(0, vl(self.f_pi, l) - vl(f, l))
        r1 = min(hf, vl(self.n_pts, l) // 2, hf)
        return (r1, (vl(self.n_pts, l) - r1))

    def _full_group_structure(self, f: int) -> tuple[int, int]:
        """Full invariants (e1, e2) of E(F_q) for curves with conductor f."""
        if f not in self.__group_cache:

            if self.is_quaternion:
                e1 = self.p**(self.n // 2) - sgn(self.t)
                e2 = e1
            else:
                e1, e2 = 1, 1
                for p, _ in factorize(self.n_pts):
                    r1, r2 = self._l_sylow(p, f)
                    e1 *= p ** r1
                    e2 *= p ** r2
            self.__group_cache[f] = (e1, e2)
        return self.__group_cache[f]

    def _torsion_subgroup(self, f: int, N: int) -> tuple[int, int]:
        """N-torsion subgroup invariants for curves with conductor f.

        Uses the cached full group structure if available (gcd path),
        otherwise computes restricted to the prime factors of N only.
        """
        if f in self.__group_cache:
            e1, e2 = self.__group_cache[f]
        else:
            # full structure not yet computed — only iterate primes of N
            e1, e2 = 1, 1
            for p, _ in factorize(N):
                r1, r2 = self._l_sylow(p, f)
                e1 *= p ** r1
                e2 *= p ** r2
        return (math.gcd(N, e1), math.gcd(N, e2))

    def _compute_torsion_record(self, N: int) -> int:
        # TODO filter

        if self.n_pts % N != 0:
            return TorsionRecord(n_curves=0, n_points=0, value=Fraction(0))

        f_list = divisors(self.f_pi)
        n_curves_total = 0
        n_pts_total = 0

        clr = Logger.SUCCESS if len(factorize(self.f_pi)) > 1 else Logger.WARNING
        if not self.ordinary:
            clr = Logger.NOTICE
        Logger.cprint(
            f"Computing torsion contribution for isogeny class with t={self.t}, N_pts={fmt_factored(self.n_pts)}, f_pi={fmt_factored(self.f_pi)}, D_K={fmt_factored(self.field.D)} ordinary={self.ordinary}",
            clr,
        )

        value_total = Fraction(0)
        for f in f_list:
            if not self.ordinary and f % self.p == 0:
                continue

            inv = self._full_group_structure(f)
            N_torsion = self._torsion_subgroup(f, N)

            n_curves = self.field.h(f)
            n_pts_exact_order = elements_of_exact_order(N, inv[0], inv[1])
            level_value = n_pts_exact_order * self.field.hw(f)

            print(
                f"f={fmt_factored(f)}, inv={fmt_invariants(inv)}, N_torsion={fmt_invariants(N_torsion)}, h={n_curves}, hw={self.field.hw(f)}, level_value={level_value}"
            )
            # for global accum
            n_curves_total += n_curves
            n_pts_total += n_curves * n_pts_exact_order
            value_total += level_value

        rec = TorsionRecord(
            n_curves=n_curves_total, n_points=n_pts_total, value=value_total
        )
        return rec


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
