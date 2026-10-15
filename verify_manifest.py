#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""results/ manifest: hashes + tally, verified or regenerated.

Format mirrors the sibling OMEGA practice so the same checker shape
applies: one line per file "SHA-256 <hex> <size> <name>", plus a
"TOTAL FILES : N" tally that must equal the row count (a manifest
that under-reports its own length is as misleading as a stale one).

Reruns rewrite results/ artifacts (they carry date fields), so any
legitimate rerun MUST regenerate the manifest afterwards:
    python verify_manifest.py --update
Failing to do so is a FAIL here, not a drift to absorb silently.

Exit codes: 0 = manifest matches scope exactly; 1 = mismatch;
            2 = results/ unreadable (never reported as pass).
"""
import hashlib
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, "results")
MANIFEST = os.path.join(RESULTS, "MANIFEST.sha256")


def scope():
    """Artifact files the manifest covers (not itself, not checkpoints)."""
    try:
        names = sorted(os.listdir(RESULTS))
    except OSError as exc:
        print("FAIL results/ unreadable: %s" % exc)
        sys.exit(2)
    return [n for n in names
            if os.path.isfile(os.path.join(RESULTS, n))
            and n != "MANIFEST.sha256"]


def digest(path):
    """sha256 + size of the LF-NORMALIZED content.

    Line-ending trap this defeats (caught by CI, not by eye): git
    checks out CRLF working copies on Windows and LF on Linux, so a
    manifest over raw bytes goes stale on every platform change even
    though every file is byte-identical in the repository. Hashing
    the normalized form makes the manifest checkout-style immune.
    """
    with open(path, "rb") as f:
        data = f.read().replace(b"\r\n", b"\n")
    h = hashlib.sha256()
    h.update(data)
    return h.hexdigest(), len(data)


def render(rows):
    lines = ["SHA-256 %s %d %s" % (h, s, n) for (n, h, s) in rows]
    lines.append("TOTAL FILES : %d" % len(rows))
    return "\n".join(lines) + "\n"


def read_manifest():
    try:
        with open(MANIFEST, "r", encoding="utf-8") as f:
            text = f.read()
    except OSError:
        return None
    rows = []
    tally = None
    for line in text.splitlines():
        if line.startswith("SHA-256 "):
            parts = line.split()
            if len(parts) != 4:
                return "malformed"
            rows.append((parts[3], parts[1], int(parts[2])))
        elif line.startswith("TOTAL FILES : "):
            tally = int(line.split(":")[1])
    if tally is None or tally != len(rows):
        return "malformed"
    return rows


def main():
    update = "--update" in sys.argv[1:]
    names = scope()
    live = []
    for n in names:
        p = os.path.join(RESULTS, n)
        h, s = digest(p)
        live.append((n, h, s))
        try:
            with open(p, "r", encoding="utf-8") as f:
                json.load(f) if p.endswith(".json") else None
        except Exception as exc:
            if p.endswith(".json"):
                print("FAIL %s is not valid JSON: %s" % (n, exc))
                return 1
        if p.endswith(".jsonl"):
            try:
                with open(p, "r", encoding="utf-8") as f:
                    for i, line in enumerate(f, 1):
                        if line.strip():
                            json.loads(line)
            except Exception as exc:
                print("FAIL %s line %d is not valid JSON: %s" % (n, i, exc))
                return 1
    if update:
        with open(MANIFEST, "w", encoding="utf-8", newline="\n") as f:
            f.write(render(live))
        print("MANIFEST regenerated: %d files" % len(live))
        return 0
    want = read_manifest()
    if want is None:
        print("FAIL MANIFEST.sha256 missing -- run with --update after a rerun")
        return 1
    if want == "malformed":
        print("FAIL MANIFEST.sha256 malformed")
        return 1
    want_map = {n: (h, s) for (n, h, s) in want}
    live_map = {n: (h, s) for (n, h, s) in live}
    if want_map != live_map:
        only_want = sorted(set(want_map) - set(live_map))
        only_live = sorted(set(live_map) - set(want_map))
        changed = sorted(n for n in set(want_map) & set(live_map)
                         if want_map[n] != live_map[n])
        if only_want:
            print("FAIL listed but absent: %s" % only_want)
        if only_live:
            print("FAIL present but unlisted: %s" % only_live)
        if changed:
            print("FAIL content drift (rerun without --update?): %s" % changed)
        return 1
    print("MANIFEST PASS: %d files, hashes match, JSON valid" % len(live))
    return 0


if __name__ == "__main__":
    sys.exit(main())
