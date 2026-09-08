import argparse
from html import escape
from pathlib import Path
from xml.etree import ElementTree

WIDTH = 200
HEIGHT = 240
DEFAULT_OUTPUT_DIR = Path("figures")

OUTLINE = "#2A8490"
CONTENT = "#438D96"
BAND = "#40969A"
PAPER = "#FFFFFF"
FOLD = "#F4FAFA"

PAPER_PATH = (
    "M 39 14 H 126 L 172 60 V 211 "
    "Q 172 222 161 222 H 39 Q 28 222 28 211 "
    "V 25 Q 28 14 39 14 Z"
)


def _svg_document(
    *,
    title: str,
    aria_label: str,
    clip_id: str,
    contents: str,
) -> str:
    return f"""<?xml version="1.0" encoding="UTF-8"?>
<svg xmlns="http://www.w3.org/2000/svg"
     width="{WIDTH}" height="{HEIGHT}" viewBox="0 0 {WIDTH} {HEIGHT}"
     role="img" aria-label="{escape(aria_label)}">
  <title>{escape(title)}</title>
  <defs>
    <clipPath id="{clip_id}">
      <path d="{PAPER_PATH}"/>
    </clipPath>
  </defs>
  <path d="{PAPER_PATH}"
        fill="{PAPER}" stroke="{OUTLINE}" stroke-width="2.6"
        stroke-linejoin="round"/>
  <path d="M 126 14 V 60 H 172 Z"
        fill="{FOLD}" stroke="{OUTLINE}" stroke-width="2.2"
        stroke-linejoin="round"/>
{contents}
</svg>
"""


def _waist_band(*, clip_id: str, label: str, font_size: int) -> str:
    return f"""  <rect x="28" y="164" width="144" height="58"
        fill="{BAND}" clip-path="url(#{clip_id})"/>
  <line x1="29.3" y1="164" x2="170.7" y2="164"
        stroke="{OUTLINE}" stroke-width="1.6"/>
  <text x="100" y="201" fill="#FFFFFF"
        font-family="Arial, Helvetica, sans-serif"
        font-size="{font_size}" font-weight="700"
        text-anchor="middle">{escape(label)}</text>"""


def generate_fasta_svg() -> str:
    clip_id = "fasta-paper-clip"
    contents = f"""  <text x="53" y="91" fill="{CONTENT}"
        font-family="Consolas, 'Courier New', monospace"
        font-size="24" font-weight="500">&gt;</text>
  <line x1="55" y1="107" x2="115" y2="107"
        stroke="{CONTENT}" stroke-width="2.8" stroke-linecap="round"/>
  <line x1="55" y1="125" x2="139" y2="125"
        stroke="{CONTENT}" stroke-width="2.8" stroke-linecap="round"/>
  <line x1="55" y1="143" x2="126" y2="143"
        stroke="{CONTENT}" stroke-width="2.8" stroke-linecap="round"/>
{_waist_band(clip_id=clip_id, label="FASTA", font_size=35)}"""
    return _svg_document(
        title="FASTA 文件图标",
        aria_label="wt.fasta file icon",
        clip_id=clip_id,
        contents=contents,
    )


def generate_json_svg() -> str:
    contents = f"""  <path d="M 77 79
        C 67 79 66 87 66 98 V 106
        C 66 116 62 121 55 123
        C 62 125 66 130 66 140 V 148
        C 66 159 67 167 77 167"
        fill="none" stroke="{CONTENT}" stroke-width="4"
        stroke-linecap="round" stroke-linejoin="round"/>
  <circle cx="100" cy="111" r="3" fill="{CONTENT}"/>
  <circle cx="100" cy="137" r="3" fill="{CONTENT}"/>
  <path d="M 123 79
        C 133 79 134 87 134 98 V 106
        C 134 116 138 121 145 123
        C 138 125 134 130 134 140 V 148
        C 134 159 133 167 123 167"
        fill="none" stroke="{CONTENT}" stroke-width="4"
        stroke-linecap="round" stroke-linejoin="round"/>"""
    return _svg_document(
        title="JSON 文件图标",
        aria_label="metadata.json file icon",
        clip_id="json-paper-clip",
        contents=contents,
    )


def generate_csv_svg() -> str:
    clip_id = "csv-paper-clip"
    contents = f"""  <rect x="53" y="82" width="70" height="64" rx="2"
        fill="{PAPER}" stroke="{CONTENT}" stroke-width="2"/>
  <rect x="53" y="82" width="70" height="14" rx="1.5"
        fill="#66A7A6"/>
  <line x1="76.3" y1="96" x2="76.3" y2="146"
        stroke="{CONTENT}" stroke-width="1.4"/>
  <line x1="99.7" y1="96" x2="99.7" y2="146"
        stroke="{CONTENT}" stroke-width="1.4"/>
  <line x1="53" y1="113" x2="123" y2="113"
        stroke="{CONTENT}" stroke-width="1.4"/>
  <line x1="53" y1="130" x2="123" y2="130"
        stroke="{CONTENT}" stroke-width="1.4"/>
{_waist_band(clip_id=clip_id, label="CSV", font_size=36)}"""
    return _svg_document(
        title="CSV 文件图标",
        aria_label="data.csv file icon",
        clip_id=clip_id,
        contents=contents,
    )


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help=f"Destination directory (default: {DEFAULT_OUTPUT_DIR})",
    )
    return parser.parse_args()


def main() -> None:
    args = _parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    icons = {
        "wt_fasta.svg": generate_fasta_svg(),
        "metadata_json.svg": generate_json_svg(),
        "data_csv.svg": generate_csv_svg(),
    }
    for filename, svg in icons.items():
        ElementTree.fromstring(svg)
        destination = args.output_dir / filename
        destination.write_text(svg, encoding="utf-8", newline="\n")
        print(f"Wrote {destination}")


if __name__ == "__main__":
    main()
