#!/usr/bin/env python3
"""Pre-compute class numbers for negative discriminants and store in SQLite.

Usage:
    python generate_classnb.py --cap 500000
    python generate_classnb.py --cap 1000000 --batch 50000
"""

import argparse
import sqlite3
import time
from pathlib import Path

from sage.libs.pari import pari

DB_PATH = Path(__file__).parent / "nt/store/classnb.db"


def init_db(conn: sqlite3.Connection) -> None:
    conn.execute("PRAGMA journal_mode=WAL")
    conn.execute("PRAGMA synchronous=NORMAL")
    conn.execute(
        "CREATE TABLE IF NOT EXISTS classnb (D INTEGER PRIMARY KEY, h INTEGER)"
    )
    conn.commit()


def is_squarefree(n: int) -> bool:
    n = abs(n)
    d = 2
    while d * d <= n:
        if n % (d * d) == 0:
            return False
        d += 1
    return True


def is_fundamental_discriminant(D: int) -> bool:
    if D >= 0:
        return False
    r = D % 4
    if r == 1:
        return is_squarefree(D)
    if r == 0:
        m = D // 4
        return m % 4 in (2, 3) and is_squarefree(m)
    return False


def generate(cap: int, batch_size: int = 20000, clear: bool = False) -> None:
    DB_PATH.parent.mkdir(parents=True, exist_ok=True)
    conn = sqlite3.connect(DB_PATH)
    if clear:
        conn.execute("DROP TABLE IF EXISTS classnb")
        conn.commit()
        print("Cleared existing table.")
    init_db(conn)

    existing = conn.execute("SELECT COUNT(*) FROM classnb").fetchone()[0]
    print(f"DB already contains {existing} entries. Starting generation up to |D|={cap}.")

    batch = []
    total = 0
    t0 = time.time()

    for D in range(-1, -(cap + 1), -1):
        if not is_fundamental_discriminant(D):
            continue
        h = int(pari(D).qfbclassno())
        batch.append((D, h))

        if len(batch) >= batch_size:
            conn.executemany("INSERT OR IGNORE INTO classnb VALUES (?,?)", batch)
            conn.commit()
            total += len(batch)
            rate = total / (time.time() - t0)
            print(f"\r  {total:>10,} stored  ({rate:,.0f}/s)  D={D}", end="", flush=True)
            batch.clear()

    if batch:
        conn.executemany("INSERT OR IGNORE INTO classnb VALUES (?,?)", batch)
        conn.commit()
        total += len(batch)

    conn.close()
    elapsed = time.time() - t0
    print(f"\nDone. {total:,} entries written in {elapsed:.1f}s  ->  {DB_PATH}")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description="Pre-compute class numbers into SQLite.")
    p.add_argument("--cap", type=int, default=1_000_000, help="Upper bound for |D|")
    p.add_argument("--batch", type=int, default=20_000, help="Batch size for commits")
    p.add_argument("--clear", action="store_true", help="Drop and recreate the table before generating")
    args = p.parse_args()
    generate(args.cap, args.batch, args.clear)
