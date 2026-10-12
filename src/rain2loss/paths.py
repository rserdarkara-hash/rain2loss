from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]   # repo root
DATA_RAW = ROOT / "data" / "raw"
DATA_PROCESSED = ROOT / "data" / "processed"
RESULTS = ROOT / "results"
LAMAH = DATA_RAW / "LamaH-CE"