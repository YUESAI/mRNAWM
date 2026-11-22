"""Download selected files from an NCBI RefSeq release directory."""

from __future__ import annotations

import argparse
import fnmatch
import html
import re
import shutil
import sys
import urllib.error
import urllib.request
from pathlib import Path


HREF_RE = re.compile(r'href=["\']([^"\']+)["\']', re.IGNORECASE)


def list_release_files(base_url: str, pattern: str) -> list[str]:
    with urllib.request.urlopen(base_url) as response:
        page = response.read().decode("utf-8", errors="replace")

    names: list[str] = []
    for match in HREF_RE.finditer(page):
        href = html.unescape(match.group(1))
        name = href.rstrip("/").split("/")[-1]
        if fnmatch.fnmatch(name, pattern):
            names.append(name)
    return sorted(set(names))


def download_file(url: str, output_path: Path, overwrite: bool = False) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    if output_path.exists() and output_path.stat().st_size > 0 and not overwrite:
        print(f"skip existing {output_path}")
        return

    tmp_path = output_path.with_suffix(output_path.suffix + ".tmp")
    print(f"download {url}")
    try:
        with urllib.request.urlopen(url) as response, tmp_path.open("wb") as handle:
            shutil.copyfileobj(response, handle)
    except urllib.error.URLError as exc:
        if tmp_path.exists():
            tmp_path.unlink()
        raise SystemExit(f"Download failed: {url}\n{exc}") from exc

    tmp_path.replace(output_path)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base-url", required=True, help="NCBI release directory URL.")
    parser.add_argument("--pattern", required=True, help="Filename glob, e.g. '*.rna.gbff.gz'.")
    parser.add_argument("--output-dir", help="Download directory.")
    parser.add_argument("--limit", type=int, default=0, help="Optional max number of files.")
    parser.add_argument("--dry-run", action="store_true", help="List matching files only.")
    parser.add_argument("--overwrite", action="store_true", help="Overwrite existing files.")
    return parser


def main() -> None:
    args = build_parser().parse_args()
    base_url = args.base_url.rstrip("/") + "/"
    names = list_release_files(base_url, args.pattern)
    if args.limit:
        names = names[: args.limit]

    if not names:
        print(f"No files matched {args.pattern!r} at {base_url}", file=sys.stderr)
        raise SystemExit(1)

    print(f"Matched {len(names)} files.")
    if args.dry_run:
        for name in names:
            print(name)
        return

    if not args.output_dir:
        raise SystemExit("--output-dir is required unless --dry-run is set.")

    output_dir = Path(args.output_dir)
    for name in names:
        download_file(base_url + name, output_dir / name, overwrite=args.overwrite)


if __name__ == "__main__":
    main()
