from __future__ import annotations
from dataclasses import dataclass
import math

from nt.elliptic_curves import CurvesRecordFq, IsogenyClass
from nt.common import legendre, factorize, divisors, euler_phi

class CCuspForm:
    def __init__(self, N: int, k: int, use_sage = False) -> None:
        self.N = N
        self.k = k
        self.level_factors = factorize(N)
        self.form = None
        self.use_sage = use_sage
        if use_sage:
            from sage.all import ZZ, CuspForms, Gamma1, ModularForms
            self.form = CuspForms(Gamma1(N), self.k) if use_sage else None

    def getInfo(self):
        if not self.use_sage:
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
        dim_S_new = S_new.dimension()
        return f"Cusp form of level {self.N} and weight {self.k}"


@dataclass
class LevelData:
    eis_term: int  # -T - sage_T
    curves_term: int  # trace computed by our implementation
    reference_val: int  # computed by sage

@dataclass
class TrFq_SGamma1Nk:
    eis_term: int  # -T - sage_T
    curves_term: int  # trace computed by our implementation
    val: int  # trace computed by our implementation
    reference_val: int  # computed by sage
    data: dict[str, int] = None

# caching here
class HeckeOperator:

    def __init__(self, cusp_form: CCuspForm) -> None:
        self.target = cusp_form
        self.data: dict[int, TrFq_SGamma1Nk] = {}  # keyed by q = p**n
        self._curves_cache: dict[int, CurvesRecordFq] = {}  # keyed by q = p**n

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

    def trace(self, p, n):
        q = p**n
        if q not in self._curves_cache:
            self._curves_cache[q] = CurvesRecordFq(p, n)
        c_req = self._curves_cache[q]
        curves_term = sum(
            I.eval_hk(self.target.k-2) * I.eval_torsion(self.target.N)
            for t, I in c_req.isogeny_classes
        )
        eis_term = self.eis_term(q)
        val = -eis_term - curves_term

        ref = None
        if self.target.form is not None:
            ref = self.target.form.hecke_operator(q).trace()

            print(self.target.getInfo())
            print(f"Trace of T_{q} on Sage gives {ref}")

        return TrFq_SGamma1Nk(
            eis_term=-eis_term,
            curves_term=curves_term,
            val=val,
            reference_val=ref,
        )
