import argparse
import math
from pathlib import Path
from typing import Sequence

Point = tuple[float, float]
Cubic = tuple[Point, Point, Point, Point]
CANVAS = 512


def _fmt(value: float) -> str:
    return f"{value:.2f}".rstrip("0").rstrip(".")


def _cubic_point(segment: Cubic, t: float) -> Point:
    p0, p1, p2, p3 = segment
    u = 1.0 - t
    return (
        u**3 * p0[0] + 3 * u**2 * t * p1[0] + 3 * u * t**2 * p2[0] + t**3 * p3[0],
        u**3 * p0[1] + 3 * u**2 * t * p1[1] + 3 * u * t**2 * p2[1] + t**3 * p3[1],
    )


def _cubic_derivative(segment: Cubic, t: float) -> Point:
    p0, p1, p2, p3 = segment
    u = 1.0 - t
    return (
        3 * u**2 * (p1[0] - p0[0]) + 6 * u * t * (p2[0] - p1[0]) + 3 * t**2 * (p3[0] - p2[0]),
        3 * u**2 * (p1[1] - p0[1]) + 6 * u * t * (p2[1] - p1[1]) + 3 * t**2 * (p3[1] - p2[1]),
    )


def _path_from_cubics(segments: Sequence[Cubic]) -> str:
    first = segments[0][0]
    commands = [f"M {_fmt(first[0])} {_fmt(first[1])}"]
    for _p0, p1, p2, p3 in segments:
        commands.append(f"C {_fmt(p1[0])} {_fmt(p1[1])}, " f"{_fmt(p2[0])} {_fmt(p2[1])}, " f"{_fmt(p3[0])} {_fmt(p3[1])}")
    return " ".join(commands)


def _path_from_points(points: Sequence[Point]) -> str:
    first, *rest = points
    commands = [f"M {_fmt(first[0])} {_fmt(first[1])}"]
    commands.extend(f"L {_fmt(x)} {_fmt(y)}" for x, y in rest)
    return " ".join(commands)


def _sample_chain(
    segments: Sequence[Cubic],
    count: int,
    *,
    start_fraction: float = 0.0,
    end_fraction: float = 1.0,
) -> list[tuple[Point, Point]]:
    dense: list[tuple[Point, Point]] = []
    samples_per_segment = 300
    for segment_index, segment in enumerate(segments):
        first_sample = 0 if segment_index == 0 else 1
        for sample_index in range(first_sample, samples_per_segment + 1):
            t = sample_index / samples_per_segment
            dense.append((_cubic_point(segment, t), _cubic_derivative(segment, t)))

    cumulative = [0.0]
    for (previous, _), (current, _) in zip(dense, dense[1:]):
        cumulative.append(cumulative[-1] + math.hypot(current[0] - previous[0], current[1] - previous[1]))

    total = cumulative[-1]
    targets = [total * (start_fraction + (end_fraction - start_fraction) * index / (count - 1)) for index in range(count)]
    result: list[tuple[Point, Point]] = []
    dense_index = 0
    for target in targets:
        while dense_index + 1 < len(cumulative) and cumulative[dense_index + 1] < target:
            dense_index += 1
        if dense_index + 1 == len(cumulative):
            result.append(dense[-1])
            continue
        span = cumulative[dense_index + 1] - cumulative[dense_index]
        ratio = 0.0 if span == 0 else (target - cumulative[dense_index]) / span
        (x1, y1), (dx1, dy1) = dense[dense_index]
        (x2, y2), (dx2, dy2) = dense[dense_index + 1]
        result.append(
            (
                (x1 + (x2 - x1) * ratio, y1 + (y2 - y1) * ratio),
                (dx1 + (dx2 - dx1) * ratio, dy1 + (dy2 - dy1) * ratio),
            )
        )
    return result


def _mix(
    first: tuple[int, int, int],
    second: tuple[int, int, int],
    amount: float,
) -> tuple[int, int, int]:
    return tuple(round(a * (1 - amount) + b * amount) for a, b in zip(first, second))


def _hex(color: tuple[int, int, int]) -> str:
    return "#" + "".join(f"{channel:02X}" for channel in color)


def _svg_document(title: str, description: str, definitions: str, body: str) -> str:
    return f"""<?xml version="1.0" encoding="UTF-8"?>
<svg xmlns="http://www.w3.org/2000/svg"
     width="{CANVAS}" height="{CANVAS}" viewBox="0 0 {CANVAS} {CANVAS}"
     fill="none" role="img" aria-labelledby="title description">
  <title id="title">{title}</title>
  <desc id="description">{description}</desc>
  <defs>
{definitions}
  </defs>
{body}
</svg>
"""


def draw_dna_svg() -> str:
    center_x = CANVAS / 2
    top_y = 36.0
    bottom_y = 476.0
    amplitude = 102.0
    turns = 2.0
    start_phase = math.pi / 2
    sample_count = 320

    def strand_point(progress: float, phase_offset: float = 0.0) -> Point:
        theta = start_phase + phase_offset + progress * turns * math.tau
        x = center_x + amplitude * math.sin(theta)
        y = top_y + (bottom_y - top_y) * progress
        return x, y

    dark_points = [strand_point(index / sample_count) for index in range(sample_count + 1)]
    light_points = [strand_point(index / sample_count, math.pi) for index in range(sample_count + 1)]
    dark_path = _path_from_points(dark_points)
    light_path = _path_from_points(light_points)

    # Each rung is an explicit complementary pair. The two coloured base units
    # meet directly at one clean junction, matching the compact reference icon.
    # Pair direction is alternated so both strands contain all four bases.
    rungs: list[str] = []
    pair_cycle = [
        ("A", "T"),
        ("G", "C"),
        ("T", "A"),
        ("C", "G"),
    ]
    base_names = {
        "A": "Adenine",
        "T": "Thymine",
        "C": "Cytosine",
        "G": "Guanine",
    }
    base_colors = {
        "A": "#1687B7",
        "T": "#91D7DF",
        "C": "#43ADC5",
        "G": "#0A4F88",
    }
    rung_count = 20
    for index in range(rung_count):
        progress = 0.035 + index * (0.93 / (rung_count - 1))
        dark_x, y = strand_point(progress)
        light_x, _ = strand_point(progress, math.pi)
        left = min(dark_x, light_x) + 8.0
        right = max(dark_x, light_x) - 8.0
        if right - left < 43.0:
            continue
        left_base, right_base = pair_cycle[index % len(pair_cycle)]
        # Lock every colour junction to the exact helix centre so all seams
        # form one clean vertical line through the icon.
        junction = center_x
        left_end = junction
        right_start = junction
        pair_code = f"{left_base}-{right_base}"
        pair_title = f"{base_names[left_base]}–{base_names[right_base]} " f"complementary base pair"
        rungs.extend(
            [
                f'    <g class="dna-base-pair" data-pair="{pair_code}" ' f'data-left-base="{left_base}" data-right-base="{right_base}" ' 'data-complementary="true">',
                f"      <title>{pair_title}</title>",
                f'      <line x1="{_fmt(left)}" y1="{_fmt(y + 1.7)}" ' f'x2="{_fmt(right)}" y2="{_fmt(y + 1.7)}" ' 'stroke="#07396F" stroke-opacity="0.48" stroke-width="12" ' 'stroke-linecap="round"/>',
                f'      <line x1="{_fmt(left)}" y1="{_fmt(y)}" ' f'x2="{_fmt(left_end)}" y2="{_fmt(y)}" ' f'stroke="{base_colors[left_base]}" ' 'stroke-width="9.5" stroke-linecap="butt"/>',
                f'      <line x1="{_fmt(right_start)}" y1="{_fmt(y)}" ' f'x2="{_fmt(right)}" y2="{_fmt(y)}" ' f'stroke="{base_colors[right_base]}" ' 'stroke-width="9.5" stroke-linecap="butt"/>',
                f'      <line x1="{_fmt(left + 3)}" y1="{_fmt(y - 1.8)}" ' f'x2="{_fmt(left_end)}" y2="{_fmt(y - 1.8)}" ' 'stroke="#E8FBFB" stroke-opacity="0.42" stroke-width="1.5" ' 'stroke-linecap="butt"/>',
                f'      <line x1="{_fmt(right_start)}" ' f'y1="{_fmt(y - 1.8)}" x2="{_fmt(right - 3)}" ' f'y2="{_fmt(y - 1.8)}" stroke="#E2FAFA" ' 'stroke-opacity="0.42" stroke-width="1.5" ' 'stroke-linecap="butt"/>',
                "    </g>",
            ]
        )

    crossings = [top_y + (bottom_y - top_y) * progress for progress in (0.125, 0.375, 0.625, 0.875)]
    clip_definitions = "\n".join(f"""    <clipPath id="dna-cross-{index}">
      <rect x="102" y="{_fmt(y - 30)}" width="308" height="60" rx="24"/>
    </clipPath>""" for index, y in enumerate(crossings))

    definitions = f"""    <linearGradient id="dna-dark-body" x1="150" y1="30"
                    x2="352" y2="480" gradientUnits="userSpaceOnUse">
      <stop offset="0" stop-color="#062D6D"/>
      <stop offset="0.34" stop-color="#0C4A91"/>
      <stop offset="0.68" stop-color="#063476"/>
      <stop offset="1" stop-color="#0A5498"/>
    </linearGradient>
    <linearGradient id="dna-light-body" x1="155" y1="30"
                    x2="357" y2="480" gradientUnits="userSpaceOnUse">
      <stop offset="0" stop-color="#168BC6"/>
      <stop offset="0.32" stop-color="#2479BA"/>
      <stop offset="0.66" stop-color="#3598CF"/>
      <stop offset="1" stop-color="#1D70AF"/>
    </linearGradient>
    <linearGradient id="dna-dark-gloss" x1="150" y1="30"
                    x2="360" y2="480" gradientUnits="userSpaceOnUse">
      <stop offset="0" stop-color="#5AA6D1" stop-opacity="0.72"/>
      <stop offset="0.55" stop-color="#2B7EB7" stop-opacity="0.5"/>
      <stop offset="1" stop-color="#78B8D8" stop-opacity="0.62"/>
    </linearGradient>
    <linearGradient id="dna-light-gloss" x1="160" y1="25"
                    x2="350" y2="485" gradientUnits="userSpaceOnUse">
      <stop offset="0" stop-color="#A4E4F0" stop-opacity="0.82"/>
      <stop offset="0.52" stop-color="#71C4DF" stop-opacity="0.58"/>
      <stop offset="1" stop-color="#8ED7E8" stop-opacity="0.7"/>
    </linearGradient>
    <filter id="dna-shadow" x="-15%" y="-10%" width="135%" height="125%">
      <feDropShadow dx="1.6" dy="2" stdDeviation="1.4"
                    flood-color="#001C46" flood-opacity="0.25"/>
    </filter>
    <path id="dna-dark-path" d="{dark_path}"/>
    <path id="dna-light-path" d="{light_path}"/>
{clip_definitions}"""

    def strand_layers(path_id: str, body: str, gloss: str) -> list[str]:
        return [
            f'    <use href="#{path_id}" stroke="#031D4F" stroke-width="32" ' 'stroke-linecap="round" stroke-linejoin="round"/>',
            f'    <use href="#{path_id}" stroke="url(#{body})" stroke-width="27" ' 'stroke-linecap="round" stroke-linejoin="round"/>',
            f'    <use href="#{path_id}" stroke="url(#{gloss})" stroke-width="5" ' 'stroke-linecap="round" stroke-linejoin="round" ' 'transform="translate(-3 -1)"/>',
        ]

    body_parts = ['  <g filter="url(#dna-shadow)">', *rungs]
    body_parts.extend(strand_layers("dna-dark-path", "dna-dark-body", "dna-dark-gloss"))
    body_parts.extend(strand_layers("dna-light-path", "dna-light-body", "dna-light-gloss"))

    # Restore alternating over/under order only inside narrow crossing windows.
    for index in (0, 2):
        body_parts.append(f'    <g clip-path="url(#dna-cross-{index})">')
        body_parts.extend(strand_layers("dna-dark-path", "dna-dark-body", "dna-dark-gloss"))
        body_parts.append("    </g>")
    for index in (1, 3):
        body_parts.append(f'    <g clip-path="url(#dna-cross-{index})">')
        body_parts.extend(strand_layers("dna-light-path", "dna-light-body", "dna-light-gloss"))
        body_parts.append("    </g>")
    body_parts.append("  </g>")

    return _svg_document(
        "Procedurally drawn DNA icon",
        ("A blue double helix generated from two mathematical sine curves, " "with directly joined alternating A-T and C-G complementary vector " "pairs, alternating strand crossings, dimensional highlights, and " "gently rounded termini."),
        definitions,
        "\n".join(body_parts),
    )


def draw_rna_svg() -> str:
    segments: list[Cubic] = [
        (
            (58.0, 392.0),
            (145.0, 426.0),
            (221.0, 373.0),
            (252.0, 284.0),
        ),
        (
            (252.0, 284.0),
            (302.0, 151.0),
            (369.0, 69.0),
            (456.0, 76.0),
        ),
    ]
    backbone = _path_from_cubics(segments)

    bases: list[str] = []
    samples = _sample_chain(
        segments,
        13,
        start_fraction=0.075,
        end_fraction=0.94,
    )
    for index, ((x, y), (dx, dy)) in enumerate(samples):
        tangent_length = math.hypot(dx, dy)
        # Clockwise normal points into the open side of this rising S curve.
        nx, ny = -dy / tangent_length, dx / tangent_length
        visible_length = 47.0 + 3.0 * math.sin(index * 1.7)
        x1, y1 = x + nx * 3.5, y + ny * 3.5
        x2, y2 = x + nx * visible_length, y + ny * visible_length
        bases.extend(
            [
                f'    <line x1="{_fmt(x1)}" y1="{_fmt(y1)}" ' f'x2="{_fmt(x2 + 0.8)}" y2="{_fmt(y2 + 1.1)}" ' 'stroke="#9F3C00" stroke-width="7.2" ' 'stroke-linecap="round"/>',
                f'    <line x1="{_fmt(x1)}" y1="{_fmt(y1)}" ' f'x2="{_fmt(x2)}" y2="{_fmt(y2)}" ' 'stroke="url(#rna-base)" stroke-width="4.8" ' 'stroke-linecap="round"/>',
            ]
        )

    definitions = """    <linearGradient id="rna-body" x1="58" y1="400"
                    x2="456" y2="66" gradientUnits="userSpaceOnUse">
      <stop offset="0" stop-color="#EF7100"/>
      <stop offset="0.3" stop-color="#E66800"/>
      <stop offset="0.72" stop-color="#CE5600"/>
      <stop offset="1" stop-color="#D65E00"/>
    </linearGradient>
    <linearGradient id="rna-face" x1="58" y1="400"
                    x2="456" y2="66" gradientUnits="userSpaceOnUse">
      <stop offset="0" stop-color="#FF8708"/>
      <stop offset="0.32" stop-color="#F47900"/>
      <stop offset="0.72" stop-color="#E26900"/>
      <stop offset="1" stop-color="#E97005"/>
    </linearGradient>
    <linearGradient id="rna-gloss" x1="70" y1="405"
                    x2="448" y2="62" gradientUnits="userSpaceOnUse">
      <stop offset="0" stop-color="#FFB44F" stop-opacity="0.82"/>
      <stop offset="0.52" stop-color="#FFA12D" stop-opacity="0.62"/>
      <stop offset="1" stop-color="#F58A18" stop-opacity="0.4"/>
    </linearGradient>
    <linearGradient id="rna-base" x1="90" y1="430"
                    x2="455" y2="75" gradientUnits="userSpaceOnUse">
      <stop offset="0" stop-color="#B94A00"/>
      <stop offset="0.45" stop-color="#CD5700"/>
      <stop offset="1" stop-color="#DB6400"/>
    </linearGradient>
    <filter id="rna-shadow" x="-10%" y="-10%" width="125%" height="125%">
      <feDropShadow dx="1.2" dy="1.6" stdDeviation="1"
                    flood-color="#7D2A00" flood-opacity="0.2"/>
    </filter>"""

    body = "\n".join(
        [
            '  <g filter="url(#rna-shadow)">',
            *bases,
            f'    <path d="{backbone}" stroke="#A74300" stroke-width="28" ' 'stroke-linecap="butt" stroke-linejoin="round"/>',
            f'    <path d="{backbone}" transform="translate(1.4 1.3)" ' 'stroke="url(#rna-body)" stroke-width="23" ' 'stroke-linecap="butt" stroke-linejoin="round"/>',
            f'    <path d="{backbone}" transform="translate(2.8 2.6)" ' 'stroke="url(#rna-face)" stroke-width="14" ' 'stroke-linecap="butt" stroke-linejoin="round"/>',
            f'    <path d="{backbone}" transform="translate(3.6 3.3)" ' 'stroke="url(#rna-gloss)" stroke-width="4" ' 'stroke-linecap="butt" stroke-linejoin="round"/>',
            "  </g>",
        ]
    )
    return _svg_document(
        "Procedurally drawn RNA icon",
        ("A dimensional orange RNA strand drawn from cubic Bezier curves " "with thirteen mathematically attached nucleotide strokes."),
        definitions,
        body,
    )


def draw_protein_svg() -> str:
    segments: list[Cubic] = [
        (
            (224.0, 87.0),
            (307.0, 52.0),
            (397.0, 83.0),
            (397.0, 165.0),
        ),
        (
            (397.0, 165.0),
            (397.0, 238.0),
            (326.0, 265.0),
            (252.0, 263.0),
        ),
        (
            (252.0, 263.0),
            (173.0, 261.0),
            (107.0, 251.0),
            (82.0, 319.0),
        ),
        (
            (82.0, 319.0),
            (58.0, 383.0),
            (108.0, 441.0),
            (174.0, 437.0),
        ),
        (
            (174.0, 437.0),
            (242.0, 433.0),
            (272.0, 359.0),
            (350.0, 369.0),
        ),
        (
            (350.0, 369.0),
            (411.0, 377.0),
            (440.0, 410.0),
            (430.0, 439.0),
        ),
    ]
    fold_path = _path_from_cubics(segments)
    residue_samples = _sample_chain(segments, 22)

    base_palette = [
        (202, 181, 235),
        (187, 158, 222),
        (173, 139, 212),
        (158, 121, 201),
    ]
    definitions = ["""    <linearGradient id="protein-link" x1="80" y1="70"
                    x2="435" y2="445" gradientUnits="userSpaceOnUse">
      <stop offset="0" stop-color="#7650AE"/>
      <stop offset="0.52" stop-color="#4D2487"/>
      <stop offset="1" stop-color="#6941A0"/>
    </linearGradient>
    <filter id="protein-shadow" x="-15%" y="-15%" width="135%" height="135%">
      <feDropShadow dx="0.8" dy="1.2" stdDeviation="1.2"
                    flood-color="#2B0D5C" flood-opacity="0.18"/>
    </filter>"""]
    shade_indices: list[int] = []
    for index in range(len(residue_samples)):
        palette_index = (index * 3 + index // 4) % len(base_palette)
        shade_indices.append(palette_index)
        center = base_palette[palette_index]
        # A small deterministic wave prevents a mechanical repeating pattern.
        variation = 0.035 * math.sin(index * 1.83)
        if variation >= 0:
            center = _mix(center, (255, 255, 255), variation)
        else:
            center = _mix(center, (46, 12, 103), -variation)
        highlight = _mix(center, (255, 255, 255), 0.46)
        shadow = _mix(center, (48, 15, 105), 0.49)
        definitions.append(f"""    <radialGradient id="residue-{index}" cx="31%" cy="24%" r="82%">
      <stop offset="0" stop-color="{_hex(highlight)}"/>
      <stop offset="0.42" stop-color="{_hex(center)}"/>
      <stop offset="1" stop-color="{_hex(shadow)}"/>
    </radialGradient>""")

    residue_parts: list[str] = []
    for index, ((x, y), _tangent) in enumerate(residue_samples):
        radius = 25.6 + 1.3 * math.sin(index * 1.37 + 0.4)
        if index in {0, len(residue_samples) - 1}:
            radius -= 0.8
        residue_parts.extend(
            [
                f'    <circle cx="{_fmt(x)}" cy="{_fmt(y)}" ' f'r="{_fmt(radius)}" fill="url(#residue-{index})" ' 'stroke="#32106F" stroke-width="4.6"/>',
                f'    <ellipse cx="{_fmt(x - radius * 0.31)}" ' f'cy="{_fmt(y - radius * 0.34)}" ' f'rx="{_fmt(radius * 0.22)}" ' f'ry="{_fmt(radius * 0.17)}" ' 'fill="#FFFFFF" opacity="0.84"/>',
            ]
        )

    body = "\n".join(
        [
            '  <g filter="url(#protein-shadow)">',
            f'    <path d="{fold_path}" stroke="url(#protein-link)" ' 'stroke-width="10" stroke-linecap="round" ' 'stroke-linejoin="round"/>',
            *residue_parts,
            "  </g>",
        ]
    )
    return _svg_document(
        "Procedurally drawn protein icon",
        ("A folded protein chain drawn from a cubic Bezier path with " "twenty-two equal-arc-length shaded amino-acid residues."),
        "\n".join(definitions),
        body,
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("figures"),
        help="destination directory (default: figures)",
    )
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    outputs = {
        "dna_icon.svg": draw_dna_svg(),
        "rna_icon.svg": draw_rna_svg(),
        "protein_icon.svg": draw_protein_svg(),
    }
    for filename, svg in outputs.items():
        destination = args.output_dir / filename
        destination.write_text(svg, encoding="utf-8", newline="\n")
        print(destination.resolve())


if __name__ == "__main__":
    main()
