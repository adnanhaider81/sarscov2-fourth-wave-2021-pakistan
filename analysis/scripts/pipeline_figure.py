#!/usr/bin/env python3

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch


STEPS = [
    {
        "key": "inputs",
        "label": "Paired-end FASTQ\nIllumina iSeq",
        "xy": (0.5, 2.1),
        "color": "#dbe7c9",
    },
    {
        "key": "qc",
        "label": "FastQC\nread assessment",
        "xy": (2.2, 2.1),
        "color": "#f4d8b0",
    },
    {
        "key": "trim",
        "label": "Trimmomatic\nadapter + Q30 trim",
        "xy": (3.9, 2.1),
        "color": "#f5b895",
    },
    {
        "key": "map",
        "label": "BWA-MEM +\nSAMtools sort/index",
        "xy": (5.6, 2.1),
        "color": "#b8d8d8",
    },
    {
        "key": "consensus",
        "label": "bcftools\nvariants + consensus",
        "xy": (7.3, 2.1),
        "color": "#9fc3d5",
    },
    {
        "key": "pangolin",
        "label": "Pangolin\nlineage assignment",
        "xy": (7.3, 0.75),
        "color": "#f8c8c8",
    },
    {
        "key": "context",
        "label": "GenBank context\naccession fetch",
        "xy": (5.6, 0.75),
        "color": "#d7c4e8",
    },
    {
        "key": "align",
        "label": "MAFFT\nwhole-genome alignment",
        "xy": (3.9, 0.75),
        "color": "#c7dce5",
    },
    {
        "key": "tree",
        "label": "IQ-TREE\nmaximum likelihood tree",
        "xy": (2.2, 0.75),
        "color": "#ffd7a8",
    },
    {
        "key": "nextstrain",
        "label": "Augur/Auspice\ninteractive build",
        "xy": (0.5, 0.75),
        "color": "#f3e7a3",
    },
]

EDGES = [
    ("inputs", "qc"),
    ("qc", "trim"),
    ("trim", "map"),
    ("map", "consensus"),
    ("consensus", "pangolin"),
    ("consensus", "context"),
    ("context", "align"),
    ("pangolin", "align"),
    ("align", "tree"),
    ("tree", "nextstrain"),
]


def add_box(ax, x, y, label, color):
    width = 1.35
    height = 0.6
    box = FancyBboxPatch(
        (x, y),
        width,
        height,
        boxstyle="round,pad=0.04,rounding_size=0.08",
        linewidth=1.6,
        edgecolor="#2a2a2a",
        facecolor=color,
    )
    ax.add_patch(box)
    ax.text(
        x + width / 2,
        y + height / 2,
        label,
        ha="center",
        va="center",
        fontsize=10,
        fontweight="bold",
        color="#202020",
    )


def add_arrow(ax, start, end):
    arrow = FancyArrowPatch(
        start,
        end,
        arrowstyle="-|>",
        mutation_scale=15,
        linewidth=1.6,
        color="#444444",
        shrinkA=10,
        shrinkB=10,
    )
    ax.add_patch(arrow)


def connect(ax, positions, source, target):
    width = 1.35
    height = 0.6
    sx, sy = positions[source]
    tx, ty = positions[target]

    if sy == ty:
        if tx > sx:
            start = (sx + width, sy + height / 2)
            end = (tx, ty + height / 2)
        else:
            start = (sx, sy + height / 2)
            end = (tx + width, ty + height / 2)
    elif ty < sy:
        start = (sx + width / 2, sy)
        end = (tx + width / 2, ty + height)
    else:
        start = (sx + width / 2, sy + height)
        end = (tx + width / 2, ty)

    add_arrow(ax, start, end)


def build_figure(output_path):
    fig, ax = plt.subplots(figsize=(11.5, 4.6), constrained_layout=True)
    fig.patch.set_facecolor("#fffdf8")
    ax.set_facecolor("#fffdf8")
    ax.set_xlim(0, 9.2)
    ax.set_ylim(0.35, 3.45)
    ax.axis("off")

    ax.text(
        0.45,
        3.18,
        "SARS-CoV-2 fourth-wave analysis pipeline",
        fontsize=18,
        fontweight="bold",
        color="#1e2a2f",
    )
    ax.text(
        0.45,
        2.94,
        "QC and consensus generation feed lineage assignment, phylogeny, and optional Nextstrain visualization.",
        fontsize=10.5,
        color="#4f5b62",
    )

    positions = {}
    for step in STEPS:
        positions[step["key"]] = step["xy"]
        add_box(ax, step["xy"][0], step["xy"][1], step["label"], step["color"])

    for source, target in EDGES:
        connect(ax, positions, source, target)

    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=220)
    plt.close(fig)


def parse_args():
    parser = argparse.ArgumentParser(description="Render a pipeline overview figure")
    parser.add_argument("--out", required=True, help="Output image path")
    return parser.parse_args()


def main():
    args = parse_args()
    build_figure(args.out)
    print(f"Saved {args.out}")


if __name__ == "__main__":
    main()
