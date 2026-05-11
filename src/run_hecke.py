from __future__ import annotations

import argparse
import random
import time

from sympy import primerange

from nt.modular_forms import CCuspForm, HeckeOperator
from utils.data import ResultData as Data
from utils.logging import Logger


def compute(p, n, N, k):
    q = p**n
    result = [Data(label="q", value=q)]
    # tr = hecke.compute_trace(p, n)
    return result


def run():
    args = parse_args()
    start_t = time.time()

    hecke = HeckeOperator(CCuspForm(args.N, args.k, args.sage))

    results = []
    if args.p == -1:
        primes = list(primerange(args.pmin, args.pmax + 1))
        if args.random:
            results.append(hecke.trace(random.choice(primes), 1))
        else:
            for p in primes:
                results.append(hecke.trace(p, 1))
    else:
        results.append(hecke.trace(args.p, args.n))

    end_t = time.time()
    q_info = f"q={args.p**args.n}" if args.p > 0 else f"prange={args.pmin}-{args.pmax}"
    Logger.header(f"Results for {q_info}  N={args.N}  k={args.k}  ({end_t - start_t:.3f}s)")
    #Logger.print_results(results)


def parse_args():
    p = argparse.ArgumentParser(description="Compute Hecke traces for S_k(Gamma_1(N)).")
    p.add_argument("-p", "--p",    type=int, default=-1, help="Field characteristic")
    p.add_argument("-n", "--n",    type=int, default=1,  help="Extension degree")
    p.add_argument("-N", "--N",    type=int, default=-1, help="Level N")
    p.add_argument("-k", "--k",    type=int, default=2,  help="Weight k")
    p.add_argument("--pmin",       type=int, default=5,  help="Lower prime bound")
    p.add_argument("--pmax",       type=int, default=100, help="Upper prime bound")
    p.add_argument("--plist",      type=int, nargs="*",  help="Explicit list of primes")
    p.add_argument("--sage",       action="store_true",  help="Compare against Sage's Hecke operator")
    p.add_argument("--random",     action="store_true",  help="Pick one random prime from prange")
    p.add_argument("--filter",     action="store_true",  help="Only process traces that can have level structure")
    p.add_argument("--flatten",    action="store_true",  help="Flatten isogeny classes into l-adic tower")
    return p.parse_args()


if __name__ == "__main__":
    run()
