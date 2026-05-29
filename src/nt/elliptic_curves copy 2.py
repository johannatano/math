from __future__ import annotations

from dataclasses import dataclass, field
from functools import cached_property, lru_cache
from fractions import Fraction
import math
from sage.libs.pari import pari
from nt.rings import ImaginaryQuadraticField, QuaternionAlgebra, _H
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
    euler_phi,
    jordan_totient,
)

from utils.logging import Logger

@dataclass
class ConductorRecord:
    """Per-conductor contribution within one isogeny class at level N."""
    f: int
    n_curves: int            # h(f), unweighted curve count
    n_pts_exact_order: int   # |{x in E(F_q)[N] : ord(x) = N}|
    full_inclusion: bool = False  # whether to include this conductor's curves in the final count (e.g. for SS, we only want to include the f=1 curves)
    value: Fraction = field(default_factory=Fraction)
    torsion_inv: tuple[int, int] = (0, 0)  # N-torsion subgroup invariants (e1, e2)
    inv: tuple[int, int] = (0, 0)  # full group invariants (e1, e2)
    n_curves_weighted: Fraction = field(default_factory=Fraction)


@dataclass
class TorsionRecord:
    """Aggregate torsion data for one isogeny class at level N."""
    n_curves: int = 0
    n_points: int = 0
    value: Fraction = field(default_factory=Fraction)
    conductor_levels: list[ConductorRecord] = field(default_factory=list)


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
            else self.field.Hf(self.f_pi)
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

    def _inert_factor(self, f: int) -> int:
        return (
            2
            if (
                not self.ordinary
                and not self.is_quaternion
                and legendre(self.field.D * f**2, self.p) < 0
            )
            else 1
        )
    def _compute_torsion_at(self, f: int, N: int, compute_full: bool = False) -> ConductorRecord:

        grp_inv = self._full_group_structure(f)
        torsion_inv = self._torsion_subgroup(f, N, compute_full)

        p_inert = self._inert_factor(f)

        _n_pts = elements_of_exact_order(N, torsion_inv[0], torsion_inv[1])
        n_curves = self.field.h_tilde(f) * p_inert
        n_curves_weighted = self.field.hw(f) * p_inert

        '''clr = Logger.ERROR if _n_pts == 0 else Logger.SUCCESS
        if vl(self.f_pi, 2) > 3 and vl(self.n_pts, 2) > 4:
            Logger.cprint(
                f"t={self.t}, f_pi={self.f_pi}, f={f}, inv={fmt_invariants(grp_inv)}, torsion_inv={fmt_invariants(torsion_inv)}, #P_{N}={_n_pts}", clr
            )'''

        has_full = torsion_inv[0] % N == 0 and torsion_inv[1] % N == 0
        clr = Logger.SUCCESS if has_full else Logger.ERROR
        '''Logger.cprint(
            f"D_K={self.field.D}, t={self.t}, f_pi={self.f_pi}, SS={not self.ordinary}, full_inclusion={has_full}, leg(DK/N)={legendre(f**2*self.field.D, N)}, torsion_inv={fmt_invariants(torsion_inv)}, group_inv={fmt_invariants(grp_inv)}, n_pts_exact_order={_n_pts}",
            clr,
        )'''

        return ConductorRecord(
            f=f,
            n_curves=n_curves if _n_pts > 0 else 0,
            n_pts_exact_order=_n_pts,
            n_curves_weighted=n_curves_weighted,
            torsion_inv=torsion_inv,
            inv=grp_inv,
            full_inclusion=has_full,
        )

    def _compute_torsion_record(self, N: int, flatten: bool = True) -> TorsionRecord:
        if N == 1:
            n_curves = self.size
            return TorsionRecord(
                n_curves=n_curves,
                n_points=1,
                value=Fraction(n_curves, self.field.u_index * 2),
                conductor_levels=[],
            )
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
        # flatten = True
        # can_flatten = flatten and not self.is_quaternion and N > 1
        # padding = 1
        '''if can_flatten:
            coprime = coprime_part(f_pi_reduced, N)
            f_pi_reduced //= coprime
            f_list = divisors(f_pi_reduced)
            if coprime > 1:
                H_coprime = self.field.Hf_inv(coprime)
        else:
            f_list = divisors(f_pi_reduced)'''

        '''conductor_levels = [
            self._compute_torsion_at(f, N) for f in f_list
        ]'''

        # has_full = any(cl.full_inclusion for cl in conductor_levels)

        hOK = self.field.O_K.h

        if self.is_quaternion:
            # BRUTE FORCE BUT ONLY ONE LEVEL HERE
            conductor_levels = [
                self._compute_torsion_at(f, N) for f in divisors(f_pi_reduced)
            ]
            return TorsionRecord(
                n_curves=sum(cl.n_curves * H_coprime * hOK for cl in conductor_levels),
                n_points=sum(
                    cl.n_pts_exact_order for cl in conductor_levels
                ),  # we only want to know exact nu pts PER curve in this level
                value=sum(
                    cl.n_pts_exact_order * cl.n_curves_weighted * H_coprime
                    for cl in conductor_levels
                ),
                conductor_levels=conductor_levels,
            )

        # if self.field.D !=-3:
        #    return TorsionRecord()

        total_sum = 1
        total_sum_ref = 1
        conductor_levels = []

        N_coprime = N

        Logger.cprint(
            f"Computing torsion record for t={self.t}, f_pi={self.f_pi}, N={N}, flatten={flatten}, coprime_part={coprime_part(self.f_pi, N)}",
            Logger.HEADLINE,
        )
        for l, a in factorize(N):

            num_non_max = 0
            num_max = 0
            # remove from coprime
            N_coprime //= l**a
            def Xi_n(i):
                return 1 if i == 0 else Fraction(l - legendre(self.field.D, l), l)
            def Psi(i):
                if a <= i <= vl(self.a_pi - 1, l):
                    return jordan_totient(l**a)
                return l ** (vl(self.f_pi, l)) * euler_phi(l**a)

            sum_at_l = 0
            sum_at_l_ref = 0
            vl_fpi = vl(self.f_pi, l)
            h = min(vl_fpi, max(0, vl(self.n_pts, l) - a))  # above max(0, vl_E-a) all zero torsion

            for i in range(0, h + 1):
                f = l ** (vl_fpi - i)
                level_sum = 0
                # level_sum = Xi_n(vl_fpi - i) * Psi(i)
                '''level_sum_ = (
                    l ** (vl(self.f_pi, l)) * euler_phi(l**a) * Xi_n(vl_fpi - i)
                )'''
                # sum_at_l_ += level_sum_
                print(
                    f"NEW: i={i}, f={f}, Xi_n={Xi_n(vl_fpi - i)}, Psi={Psi(i)}, level_sum_={level_sum}"
                )
                # sum_at_l += level_sum

            for i in range(0, vl_fpi + 1):
                f = l ** (vl_fpi - i)
                if (
                    i > vl(self.n_pts, l) - a
                ):  # everything above this conductor level has zero l**a torsion, since vl(e1) < a
                    print(
                        f"i={i}, f={f}, EXIT ZERO"
                    )
                    break
                if i == a:
                    level_sum = jordan_totient(l**a) * self.field.Hf_tilde(f)
                    sum_at_l_ref += level_sum
                    num_max = self.field.Hf_tilde(f)
                    break
                if i == vl(self.a_pi - 1, l):
                    # we reached max torsion allowed
                    level_sum = l**i * euler_phi(l**a) * self.field.Hf_tilde(f)
                    sum_at_l_ref += level_sum
                    num_non_max += self.field.Hf_tilde(f)
                    print(
                        f"i={i}, f={f}, adding CAP level={level_sum}, num_non_max={num_non_max}"
                    )
                    break
                num_non_max += 1
                xi = Fraction(l - legendre(self.field.D, l), l) if f > 1 else Fraction(1)
                #level_sum = l ** (vl(self.f_pi, l)) * euler_phi(l**a) * xi

                level_sum = l**i * euler_phi(l**a) * self.field.h_tilde(f)
                sum_at_l_ref += level_sum
                print(
                    f"i={i}, f={f}, adding ONE level_sum_={level_sum}"
                )

            total_sum_ref *= sum_at_l_ref
            total_sum *= sum_at_l_ref
            # level_sum = l**i * euler_phi(l**a) * self.field.h_tilde(f)
            # sum_at_l += level_sum'''
        print(f"fpi={self.f_pi}, N={N}, N_coprime={N_coprime}, total_sum={total_sum}")
        H_coprime_ = self.field.Hf_tilde(N_coprime)
        final_value = total_sum_ref * Fraction(hOK * H_coprime_, self.field.u_index * 2)

        return TorsionRecord(
            n_curves=0,  # sum(cl.n_curves for cl in conductor_levels) * H_coprime_ * hOK,
            n_points=0,  # sum(
            # cl.n_pts_exact_order for cl in conductor_levels
            # ),  # we only want to know exact nu pts PER curve in this level
            value=final_value,
            conductor_levels=[],  # conductor_levels,
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
            ts.append((0, None))
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
        return Fraction(0) + sum(I.size for t, I in self.isogeny_classes)

    @cached_property
    def isogeny_classes(self):
        return self.__isogeny_classes.items()

    def get_isogeny_class(self, t: int):
        return self.__isogeny_classes.get(t)
