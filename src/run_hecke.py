from __future__ import annotations

import argparse
import random
import time

from sympy import primerange

from nt.modular_forms import CuspForm, HeckeOperator
from nt.common import fmt_magnitude, equiv
from utils.data import ResultData as Data
from utils.logging import Logger
from math import gcd

def format_result_data(data: TrFq_SGamma1Nk, level: int) -> list[Data]:
    return [
        Data("p", data.p),
        Data("q", data.q),
        Data("(q,N)", gcd(data.q, level)),
        Data("q mod N", equiv(data.q, level)),
        Data("Cusps", data.eis_term),
        Data("Y1(N)", data.curves_term),
        Data("Computed", data.val),
        Data("Ref", data.reference_val),
        Data("Error", data.error, fmt="factors"),  # , fmt="factors"
        Data("#E", data.num_curves),
        Data("#E(SS)", data.num_ss_curves),
        Data("#P", data.num_points),
    ]

def run():
    args = parse_args()
    start_t = time.time()
    hecke_op = HeckeOperator(CuspForm(args.N, args.k+2, args.sage))

    results = []
    if args.p == -1:
        primes = list(primerange(args.pmin, args.pmax + 1))
        if args.random:
            results.append(format_result_data(hecke_op.trace(random.choice(primes), 1), args.N))
        else:
            for p in primes:
               results.append(format_result_data(hecke_op.trace(p, args.n), args.N))
    else:
        results.append(format_result_data(hecke_op.trace(args.p, args.n), args.N))

    end_t = time.time()
    
    q = args.p**args.n if args.p > 0 else f"{args.pmin}-{args.pmax}"
    q_info = f"q={q} magnitude={fmt_magnitude(q)}" if args.p > 0 else f"prange={args.pmin}-{args.pmax}"
    Logger.header(f"Results for {q_info}  N={args.N}  k={args.k}  ({end_t - start_t:.3f}s)")
    Logger.print_results(results)


def parse_args():
    p = argparse.ArgumentParser(description="Compute Hecke traces for S_k(Gamma_1(N)).")
    p.add_argument("-p", "--p",    type=int, default=-1, help="Field characteristic")
    p.add_argument("-n", "--n",    type=int, default=1,  help="Extension degree")
    p.add_argument("-N", "--N",    type=int, default=11, help="Level N")
    p.add_argument("-k", "--k",    type=int, default=1,  help="Weight k")
    p.add_argument("--pmin",       type=int, default=2,  help="Lower prime bound")
    p.add_argument("--pmax",       type=int, default=100, help="Upper prime bound")
    p.add_argument("--plist",      type=int, nargs="*",  help="Explicit list of primes")
    p.add_argument("--sage",       action="store_true",  help="Compare against Sage's Hecke operator")
    p.add_argument("--random",     action="store_true",  help="Pick one random prime from prange")
    p.add_argument("--filter",     action="store_true",  help="Only process traces that can have level structure")
    p.add_argument("--flatten",    action="store_true",  help="Flatten isogeny classes into l-adic tower")
    return p.parse_args()


if __name__ == "__main__":
    run()
