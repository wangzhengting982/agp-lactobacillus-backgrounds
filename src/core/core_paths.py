from pathlib import Path
ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data/archive/01_核心统计"
OUT = ROOT / "runs/core"
OUT.mkdir(parents=True, exist_ok=True)
PRIMARY = ROOT / 'runs/preparation/project/03_分析与结果/association_models'
TAXONOMY_INPUT = OUT / 'taxonomy_inputs'
