"""Render both current and submitted CMYK PDFs at 200 dpi and compare pixels.

Exit zero only when all three single-page PDFs have exactly identical RGB pixels.
This verifies artwork rendering; statistical reproduction is checked separately.
Neither pre-existing preview images nor PDF file hashes substitute for rendering.
"""
from pathlib import Path
import argparse
import hashlib
import importlib.metadata
import json
import sys
from datetime import datetime, timezone

import numpy as np
from PIL import ImageChops
import pypdfium2 as pdfium

ROOT = Path(__file__).resolve().parents[1]
NAMES = [
    '图1_分析人群与乳酸菌检出.pdf',
    '图2_属间比较与饮食调整.pdf',
    '图3_背景水平与检出及丰度.pdf',
]
DPI = 200


def resolve(value):
    value = Path(value)
    return value.resolve() if value.is_absolute() else (ROOT / value).resolve()


def relative(value):
    try:
        return value.relative_to(ROOT).as_posix()
    except ValueError:
        return str(value)


def sha256(value):
    return hashlib.sha256(value.read_bytes()).hexdigest()


def render(value):
    doc = pdfium.PdfDocument(str(value))
    try:
        if len(doc) != 1:
            raise ValueError(f'Expected one page, found {len(doc)}')
        page = doc[0]
        try:
            bitmap = page.render(scale=DPI / 72)
            try:
                return bitmap.to_pil().convert('RGB').copy()
            finally:
                bitmap.close()
        finally:
            page.close()
    finally:
        doc.close()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--current-dir', type=Path, default=Path('runs/figures/vector'))
    parser.add_argument('--reference-dir', type=Path,
                        default=Path('data/archive/11_投稿格式修订_20261006/图件复绘/vector'))
    parser.add_argument('--report', type=Path, default=Path('reports/figures_pixel_comparison.json'))
    args = parser.parse_args()
    current, reference, destination = map(resolve, [args.current_dir, args.reference_dir, args.report])
    checks = []
    for name in NAMES:
        new, old = current / name, reference / name
        result = {'figure': Path(name).stem, 'current_pdf': relative(new),
                  'submitted_pdf': relative(old), 'passed': False}
        try:
            result.update({'current_pdf_sha256': sha256(new),
                           'submitted_pdf_sha256': sha256(old)})
            a, b = render(new), render(old)
            try:
                result.update({'current_size_pixels': list(a.size),
                               'submitted_size_pixels': list(b.size)})
                if a.size != b.size:
                    result['error'] = 'Rendered pixel dimensions differ'
                else:
                    difference = ImageChops.difference(a, b)
                    pixels = np.asarray(difference)
                    changed = int(np.any(pixels != 0, axis=2).sum())
                    result.update({'differing_pixels': changed,
                                   'maximum_channel_difference': int(pixels.max()),
                                   'difference_bbox': difference.getbbox(),
                                   'pixel_identical_to_submitted_vector_render': changed == 0,
                                   'passed': changed == 0})
                    difference.close()
            finally:
                a.close()
                b.close()
        except Exception as error:
            result['error'] = f'{type(error).__name__}: {error}'
        checks.append(result)
    passed = len(checks) == 3 and all(item['passed'] for item in checks)
    report = {'status': 'PASS' if passed else 'FAIL', 'all_passed': passed,
              'checked_at_utc': datetime.now(timezone.utc).isoformat(),
              'dpi': DPI, 'render_scale': DPI / 72, 'comparison_mode': 'RGB exact pixel equality',
              'tolerance': 0, 'figures_checked': len(checks),
              'pypdfium2_version': importlib.metadata.version('pypdfium2'),
              'pillow_version': importlib.metadata.version('Pillow'),
              'checks': checks,
              'scope': 'Both PDF sets are freshly rendered in the same process. No model fitting or independent visual inspection is performed by this check.'}
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding='utf-8')
    print(json.dumps({'status': report['status'], 'figures_checked': len(checks),
                      'report': relative(destination)}, ensure_ascii=False))
    return 0 if passed else 1


if __name__ == '__main__':
    if hasattr(sys.stdout, 'reconfigure'):
        sys.stdout.reconfigure(encoding='utf-8')
    raise SystemExit(main())
