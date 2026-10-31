"""Download the two external inputs. Run once, from the repo root.

    python scripts/fetch_inputs.py

Writes data/1n8z.pdb and data/her2_ectodomain.fasta (P04626, Thr23-Thr652).
Both are gitignored - everyone fetches their own copy.
"""
import urllib.request
from pathlib import Path

DATA = Path(__file__).resolve().parent.parent / "data"

def get(url: str, dest: Path) -> str:
    print(f"  {url}")
    text = urllib.request.urlopen(url, timeout=60).read().decode()
    dest.write_text(text)
    return text

print("1N8Z (trastuzumab Fab + HER2, 2.52 A)")
get("https://files.rcsb.org/download/1N8Z.pdb", DATA / "1n8z.pdb")

print("HER2 P04626")
raw = get("https://rest.uniprot.org/uniprotkb/P04626.fasta", DATA / "_p04626_full.fasta")
full = "".join(l for l in raw.splitlines() if not l.startswith(">"))
ecto = full[22:652]          # Thr23 - Thr652, 1-based inclusive
assert ecto[0] == "T" and ecto[-1] == "T", f"boundaries off: {ecto[0]}..{ecto[-1]}"
(DATA / "her2_ectodomain.fasta").write_text(
    f">HER2_ectodomain_P04626_T23-T652 {len(ecto)}aa\n{ecto}\n"
)
print(f"  ectodomain {len(ecto)} aa, {ecto[:8]}...{ecto[-8:]}")
print("\ndone")
