#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Endpoint-convention test: what happens when mu hits an eigenvalue.

CONVENTION UNDER TEST (see NOTES.md C.8): N(mu) = n_neg counts
eigenvalues STRICTLY below mu. An eigenvalue exactly AT mu must
never be assigned to either side: the LDL^T pivot goes exactly zero
and the count comes back UNDETERMINED (fail-closed), never guessed.

Matrix: diag(1, 2) -- eigenvalues exactly 1 and 2, no arithmetic
doubt possible. Requirements:
  mu = 0.999999 -> n_neg == 0, determined (nothing strictly below 1)
  mu = 1.0      -> UNDETERMINED (eigenvalue on the shift; unassigned)
  mu = 1.5      -> n_neg == 1, determined (strictly-below semantics)
  mu = 2.0      -> UNDETERMINED (same, second pivot)
  mu = 2.000001 -> n_neg == 2, determined

Needs ../GUINAND_WEIL/source_arb_ldlt_certify.py (same requirement as
every other test here). Without it: loud SKIP, exit 0 -- never a
silent pass and never a failure for a missing optional tree.
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

TREE = os.path.normpath(os.path.join(HERE, "..", "GUINAND_WEIL",
                                      "source_arb_ldlt_certify.py"))

if not os.path.exists(TREE):
    print("SKIP test_endpoint_convention: ../GUINAND_WEIL tree absent "
          "(local-only test; CI runs py_compile + manifest instead)")
    sys.exit(0)

from sturm_count import count_shifted, load_source  # noqa: E402

from flint import arb_mat  # noqa: E402

os.environ.setdefault("OMP_NUM_THREADS", "1")


def main():
    src = load_source()
    A = arb_mat([[1, 0], [0, 2]])
    failures = []

    def check(mu, want_neg=None, want_undet=False):
        rec = count_shifted(src, A, 2, mu)
        undet = rec["undetermined"] is not None
        if want_undet:
            ok = undet
            detail = "undetermined=%r" % rec["undetermined"]
        else:
            ok = (not undet) and rec["n_neg"] == want_neg
            detail = "n_neg=%r undetermined=%r" % (rec["n_neg"],
                                                  rec["undetermined"])
        print("  mu=%-10s -> %s %s" % (mu, detail,
                                       "ok" if ok else "FAIL"), flush=True)
        if not ok:
            failures.append(mu)

    check(0.999999, want_neg=0)
    check(1.0, want_undet=True)
    check(1.5, want_neg=1)
    check(2.0, want_undet=True)
    check(2.000001, want_neg=2)

    if failures:
        print("ENDPOINT CONVENTION: FAIL at mu=%r" % failures, flush=True)
        return 1
    print("ENDPOINT CONVENTION: PASS (exact-hit never assigned; "
          "strict-below exact)", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
