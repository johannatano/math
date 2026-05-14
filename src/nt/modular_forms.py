from __future__ import annotations
from dataclasses import dataclass
import math
from nt.elliptic_curves import CurvesRecordFq, IsogenyClass
from nt.common import legendre, factorize, divisors, euler_phi
from utils.logging import Logger


class CuspForm:
    def __init__(self, N: int, k: int, use_sage=False) -> None:
        self.N = N
        self.k = k
        self.level_factors = factorize(N)
        self.form = None
        self.use_sage = use_sage
        if use_sage:
            from sage.all import ZZ, CuspForms, Gamma1, ModularForms

            self.form = CuspForms(Gamma1(N), self.k) if use_sage else None

    def getInfo(self):
        return f"Cusp form"
        """if not self.use_sage:
            return f"No Sage form for level {self.N} and weight {self.k}"
        from sage.all import ZZ, CuspForms, Gamma1, ModularForms
        M = ModularForms(Gamma1(self.N), self.k)
        E = M.eisenstein_subspace()
        # EISEN IS ALWAYS #CUSPS
        S = M.cuspidal_subspace()
        dim_M = M.dimension()
        dim_E = E.dimension()
        dim_S = S.dimension()
        S_old = S.old_submodule()
        S_new = S.new_submodule()
        dim_S_old = S_old.dimension()
        dim_S_new = S_new.dimension()"""
        # return f"Cusp form of level {self.N} and weight {self.k}"


@dataclass
class LevelData:
    eis_term: int  # -T - sage_T
    curves_term: int  # trace computed by our implementation
    reference_val: int  # computed by sage


@dataclass
class TrFq_SGamma1Nk:
    eis_term: int = 0
    curves_term: int = 0
    val: int = 0
    reference_val: int = 0
    error: int = 0
    num_curves: int = 0
    num_ss_curves: int = 0
    num_points: int = 0
    q: int = 0
    traces: list[int] = None


# caching here
class HeckeOperator:

    def __init__(self, cusp_form: CuspForm) -> None:
        self.target = cusp_form
        self.data: dict[int, TrFq_SGamma1Nk] = {}
        self._curves_cache: dict[int, CurvesRecordFq] = {}  # survives target changes

    def update_target(self, cusp_form: CuspForm) -> None:
        """Change level/weight without discarding the curves cache."""
        self.target = cusp_form
        self.data.clear()

    def eis_term(self, q):
        split = 0
        non_split = 0
        k = self.target.k - 2
        for d in divisors(self.target.N):
            Nd = self.target.N // d
            # split cusps: need (N/d) | (q-1)
            if (q - 1) % Nd == 0:
                split += euler_phi(d) * euler_phi(Nd)
            # non-split cusps: need d | 2 AND (N/d) | (q+1)
            if d in (1, 2) and (q + 1) % Nd == 0:
                non_split += euler_phi(d) * euler_phi(Nd)
        return (split * 1 + non_split * ((-1) ** self.target.k)) // 2

    def get_isogeny_class(self, q: int, t: int):
        if q not in self._curves_cache:
            return None
        return self._curves_cache[q].get_isogeny_class(t)

    def trace(self, p, n):
        q = p**n
        if q > 10**10:
            print(
                f"Too large q={q}, skipping curve computations and returning 0 for trace."
            )
            return TrFq_SGamma1Nk(
                eis_term=0,
                curves_term=0,
                val=0,
                reference_val=0,
                error=0,
                num_curves=0,
                num_ss_curves=0,
                q=q,
            )

        if q not in self._curves_cache:
            self._curves_cache[q] = CurvesRecordFq(p, n)

        c_req = self._curves_cache[q]

        N, k = self.target.N, self.target.k
        curves_term = 0
        num_curves = 0
        num_points = 0
        num_ss_curves = 0
        traces = []
        for _, I in c_req.isogeny_classes:
            if I.n_pts % N != 0:
                continue
            tor = I.get_torsion(N)
            curves_term += I.eval_hk(k - 2) * tor.value
            num_curves += tor.n_curves
            num_points += tor.n_points
            num_ss_curves += tor.n_curves if not I.ordinary else 0
            traces.append(_)

        eis_term = self.eis_term(q)
        val = -eis_term - curves_term

        ref = 0
        if self.target.form is not None:
            ref = int(self.target.form.hecke_operator(q).trace())

        return TrFq_SGamma1Nk(
            eis_term=-eis_term,
            curves_term=curves_term,
            val=val,
            reference_val=ref,
            error=ref - val,
            num_curves=num_curves,
            num_ss_curves=num_ss_curves,
            num_points=num_points,
            traces=traces,
            q=q,
        )
