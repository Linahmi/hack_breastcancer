"""Test run_cofold.py's pure logic with modal stubbed out.

Usage: python3 -I test_runner.py <path to src/aida>
"""
import sys
import types
from pathlib import Path

AIDA = Path(sys.argv[1]).resolve()

# --- stub modal_boltz2 so run_cofold imports without touching Modal ---
stub = types.ModuleType("modal_boltz2")


class _App:
    def local_entrypoint(self):
        def deco(fn):
            return fn
        return deco

    def function(self, **kw):
        def deco(fn):
            return fn
        return deco


stub.app = _App()
stub.boltz = lambda *a, **k: None
sys.modules["modal_boltz2"] = stub

sys.path.insert(0, str(AIDA))
import run_cofold as R  # noqa: E402

PARENT = (
    "EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFT"
    "ISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSSGGGGSGGGGSGGGGSDIQ"
    "MTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFT"
    "LTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK"
)

fails = []


def check(label, fn, want_ok, expect_sub=None):
    try:
        fn()
        ok = True
        msg = ""
    except Exception as e:
        ok = False
        msg = f"{type(e).__name__}: {e}"
    if ok != want_ok:
        fails.append(f"{label}: wanted {'ok' if want_ok else 'raise'}, got {msg or 'ok'}")
    elif expect_sub and expect_sub.lower() not in msg.lower():
        fails.append(f"{label}: raised but message lacked {expect_sub!r}: {msg}")
    print(f"  {'ok  ' if ok == want_ok else 'FAIL'} {label}"
          + (f"  <- {msg}" if msg else ""))


print("=== construct loads ===")
her2, b_ss = R.load_her2()
print(f"  HER2 {len(her2)} aa, {len(b_ss)} disulfide pairs, local {b_ss[0]}..{b_ss[-1]}")
assert len(her2) == 142, her2
assert len(b_ss) == 10

print("\n=== run name rules ===")
check("plain name", lambda: R.check_run_name("design_001"), True)
check("rejects 'pae' substring", lambda: R.check_run_name("paex1"), False, "pae")
check("rejects 'PAE' uppercase", lambda: R.check_run_name("myPAEtest"), False, "pae")
check("rejects slash", lambda: R.check_run_name("a/b"), False, "unsafe")
check("rejects space", lambda: R.check_run_name("a b"), False, "unsafe")
check("rejects empty", lambda: R.check_run_name(""), False, "empty")
check("allows seed suffix", lambda: R.check_run_name("parent_seed3"), True)

print("\n=== row validation ===")
check("parent valid", lambda: R.validate_row("parent", PARENT), True)
check("lowercase coerced", lambda: R.validate_row("p2", PARENT.lower()), True)
check("too long", lambda: R.validate_row("p3", PARENT + "A" * 20), False, "260")
check("bad residue X", lambda: R.validate_row("p4", PARENT[:-1] + "X"), False, "non-standard")
check("extra cysteine", lambda: R.validate_row("p5", PARENT[:-1] + "C"), False, "expected exactly 4")
check("lost cysteine", lambda: R.validate_row("p6", PARENT.replace("C", "S", 1)), False, "expected exactly 4")
check("linker deleted", lambda: R.validate_row("p7", PARENT.replace(R.LINKER, "")), False, "linker")
check("empty seq", lambda: R.validate_row("p8", ""), False, "empty")

print("\n=== disulfide pairing follows insertions ===")
# insert 6 aa into CDR-H3 region (after pos 100), pushing VL cysteines +6
shifted = PARENT[:100] + "GSGSGS" + PARENT[100:]
seq, a_ss = R.validate_row("ins6", shifted)
print(f"  len {len(seq)}  pairs {a_ss}")
assert a_ss == [(22, 96), (164, 229)], a_ss
print("  ok   pairs tracked the insertion (158->164, 223->229)")

print("\n=== YAML emission ===")
seq, a_ss = R.validate_row("parent", PARENT)
y = R.build_yaml(seq, a_ss, her2, b_ss)
assert y.count("- bond:") == 12, y.count("- bond:")
assert y.index("id: A") < y.index("id: B"), "chain order wrong"
assert "atom1: [A, 22, SG]" in y
assert "atom1: [B, 1, SG]" in y and "atom2: [B, 132, SG]" in y
assert "pae" not in y
print(f"  ok   12 bonds, chain A before B, {len(y)} bytes")

print("\n=== output path contract ===")
p = R.expected_npz(Path("outputs"), "design_001")
print(f"  {p}")
assert str(p) == ("outputs/design_001/boltz_results_design_001/predictions/"
                  "design_001/pae_design_001_model_0.npz"), p
print("  ok   matches the layout verified on a real boltz 2.2.1 run")

print("\n=== BASE_PARAMS ===")
assert "--write_full_pae" in R.BASE_PARAMS
print(f"  ok   --write_full_pae present: {R.BASE_PARAMS}")

print()
if fails:
    print(f"{len(fails)} FAILURE(S):")
    for f in fails:
        print("  " + f)
    raise SystemExit(1)
print("all checks passed")
