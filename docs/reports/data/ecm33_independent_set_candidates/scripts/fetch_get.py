# ruff: noqa
# Scoping script kept exactly as run on 2026-10-23 (stdlib only).
import sys, urllib.request

u = sys.argv[1]
r = urllib.request.urlopen(
    urllib.request.Request(u, headers={"User-Agent": "adhesionPred-scoping/1.0", "Accept": "*/*"}),
    timeout=60,
)
sys.stdout.write(r.read().decode())
