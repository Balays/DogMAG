#!/usr/bin/env python3
"""Create the publication figure for paired RefSeq versus DogMAG recruitment."""

from __future__ import annotations

import csv
import hashlib
import html
import statistics
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
INPUT = ROOT / "waltham_mapping_pilot/minitax_DogMAG/results/final_balanced46_comparison/detailed_per_sample_metrics.tsv"
OUTPUT = ROOT / "figures/dogmag_final_20260727/figure6_refseq_vs_dogmag_mapping_20260729.svg"

WIDTH, HEIGHT = 1600, 1030
REFSEQ = "#3B6FB6"
DOGMAG = "#D95F02"
GRID = "#D9DDE3"
TEXT = "#20242A"
MUTED = "#66707A"
PAIR = "#A8AFB8"


def read_rows():
    with INPUT.open(newline="", encoding="utf-8-sig") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    for row in rows:
        for key in list(row):
            if key.endswith("fraction_of_primary") or key.startswith("delta_"):
                row[key] = float(row[key])
    return rows


def quantile(values, q):
    ordered = sorted(values)
    pos = (len(ordered) - 1) * q
    lo = int(pos)
    hi = min(lo + 1, len(ordered) - 1)
    frac = pos - lo
    return ordered[lo] * (1 - frac) + ordered[hi] * frac


def esc(value):
    return html.escape(str(value))


def jitter(sample):
    digest = hashlib.md5(sample.encode("utf-8")).hexdigest()
    return (int(digest[:8], 16) / 0xFFFFFFFF - 0.5) * 42


def main():
    rows = read_rows()
    panels = [
        ("AllKennel", "Primary reads mapped", "mapped_fraction_of_primary", "6.60e-5", "A"),
        ("AllKennel", "Primary reads assigned", "assigned_fraction_of_primary", "5.72e-6", "B"),
        ("Waltham", "Primary reads mapped", "mapped_fraction_of_primary", "2.38e-7", "C"),
        ("Waltham", "Primary reads assigned", "assigned_fraction_of_primary", "2.38e-7", "D"),
    ]
    lefts = [150, 850]
    tops = [90, 555]
    pw, ph = 590, 380
    plot_top_offset, plot_bottom_offset = 72, 315
    x_ref_offset, x_dog_offset = 205, 415

    svg = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{WIDTH}" height="{HEIGHT}" viewBox="0 0 {WIDTH} {HEIGHT}">',
        '<rect width="100%" height="100%" fill="#FFFFFF"/>',
        '<style>text{font-family:Arial,Helvetica,sans-serif;fill:#20242A}.panel{font-size:24px;font-weight:700}.metric{font-size:22px;font-weight:700}.axis{font-size:17px}.small{font-size:16px}.median{font-size:17px;font-weight:700}.pval{font-size:17px;font-style:italic}</style>',
        f'<circle cx="610" cy="42" r="8" fill="{REFSEQ}"/><text x="627" y="48" class="axis">Bacterial RefSeq</text>',
        f'<circle cx="795" cy="42" r="8" fill="{DOGMAG}"/><text x="812" y="48" class="axis">DogMAG 95% ANI</text>',
        f'<line x1="1015" y1="42" x2="1050" y2="42" stroke="{PAIR}" stroke-width="2"/><text x="1062" y="48" class="axis">Paired library</text>',
    ]

    for idx, (cohort, metric_label, metric, pvalue, letter) in enumerate(panels):
        col, row_idx = idx % 2, idx // 2
        left, top = lefts[col], tops[row_idx]
        plot_top = top + plot_top_offset
        plot_bottom = top + plot_bottom_offset
        x_ref, x_dog = left + x_ref_offset, left + x_dog_offset

        members = [r for r in rows if r["cohort"] == cohort]
        ref_values = [r[f"refseq_{metric}"] for r in members]
        dog_values = [r[f"dogmag_{metric}"] for r in members]
        ref_med, dog_med = statistics.median(ref_values), statistics.median(dog_values)
        ref_q1, ref_q3 = quantile(ref_values, 0.25), quantile(ref_values, 0.75)
        dog_q1, dog_q3 = quantile(dog_values, 0.25), quantile(dog_values, 0.75)

        def y(value):
            return plot_bottom - value * (plot_bottom - plot_top)

        svg += [
            f'<text x="{left}" y="{top+20}" class="panel">{letter}</text>',
            f'<text x="{left+38}" y="{top+20}" class="metric">{esc(metric_label)}</text>',
            f'<text x="{left+38}" y="{top+48}" class="small" fill="{MUTED}">{"Mixed-kennel source cohort" if cohort == "AllKennel" else "Independent Waltham cohort"} (n=23)</text>',
            f'<text x="{left+pw-5}" y="{top+25}" text-anchor="end" class="pval">exact sign test: P = {pvalue}</text>',
        ]

        for tick in range(0, 101, 20):
            yy = y(tick / 100)
            svg.append(f'<line x1="{left+95}" y1="{yy:.1f}" x2="{left+pw-15}" y2="{yy:.1f}" stroke="{GRID}" stroke-width="1"/>')
            svg.append(f'<text x="{left+82}" y="{yy+6:.1f}" text-anchor="end" class="axis">{tick}</text>')
        svg.append(f'<line x1="{left+95}" y1="{plot_top}" x2="{left+95}" y2="{plot_bottom}" stroke="{TEXT}" stroke-width="1.5"/>')
        if col == 0:
            svg.append(f'<text transform="translate({left+25},{(plot_top+plot_bottom)/2}) rotate(-90)" text-anchor="middle" class="axis">Primary reads (%)</text>')

        for rec in members:
            j = jitter(rec["sample"])
            rv, dv = rec[f"refseq_{metric}"], rec[f"dogmag_{metric}"]
            svg.append(f'<line x1="{x_ref+j:.1f}" y1="{y(rv):.1f}" x2="{x_dog+j:.1f}" y2="{y(dv):.1f}" stroke="{PAIR}" stroke-width="1.6" opacity="0.62"/>')
            svg.append(f'<circle cx="{x_ref+j:.1f}" cy="{y(rv):.1f}" r="4.4" fill="{REFSEQ}" opacity="0.76"/>')
            svg.append(f'<circle cx="{x_dog+j:.1f}" cy="{y(dv):.1f}" r="4.4" fill="{DOGMAG}" opacity="0.76"/>')

        for xpos, med, q1, q3, color in (
            (x_ref, ref_med, ref_q1, ref_q3, REFSEQ),
            (x_dog, dog_med, dog_q1, dog_q3, DOGMAG),
        ):
            svg.append(f'<line x1="{xpos}" y1="{y(q1):.1f}" x2="{xpos}" y2="{y(q3):.1f}" stroke="{TEXT}" stroke-width="8" stroke-linecap="round"/>')
            yy = y(med)
            pts = f'{xpos},{yy-10:.1f} {xpos+10},{yy:.1f} {xpos},{yy+10:.1f} {xpos-10},{yy:.1f}'
            svg.append(f'<polygon points="{pts}" fill="{color}" stroke="#FFFFFF" stroke-width="2"/>')
            svg.append(f'<text x="{xpos}" y="{plot_bottom+28}" text-anchor="middle" class="median">{100*med:.2f}%</text>')

        svg += [
            f'<text x="{x_ref}" y="{plot_bottom+57}" text-anchor="middle" class="axis">RefSeq</text>',
            f'<text x="{x_dog}" y="{plot_bottom+57}" text-anchor="middle" class="axis">DogMAG</text>',
        ]

    svg.append('</svg>')
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text("\n".join(svg), encoding="utf-8")
    print(OUTPUT)


if __name__ == "__main__":
    main()
