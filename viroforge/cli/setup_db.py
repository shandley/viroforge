#!/usr/bin/env python3
"""
ViroForge database setup command.

Downloads required data files and builds the local database.

Usage:
    viroforge setup-db                    # Full setup (download + build + taxonomy)
    viroforge setup-db --download-only    # Only download raw data
    viroforge setup-db --skip-taxonomy    # Build database without ICTV taxonomy

Author: ViroForge Development Team
"""

import logging
import sys
import urllib.request
import shutil
from pathlib import Path

logger = logging.getLogger(__name__)

# ICTV VMR download URL (MSL40 - matches parse_ictv_taxonomy.py sheet name)
ICTV_VMR_URL = "https://ictv.global/sites/default/files/VMR/VMR_MSL40.v2.20260223.xlsx"
ICTV_VMR_FILENAME = "VMR_current.xlsx"


def download_with_progress(url: str, dest: Path, label: str = "Downloading") -> None:
    """Download a file with a progress indicator."""
    print(f"{label}: {url}")
    print(f"  -> {dest}")

    try:
        req = urllib.request.Request(url, headers={"User-Agent": "ViroForge/0.12.0"})
        with urllib.request.urlopen(req) as response:
            total = response.getheader("Content-Length")
            total = int(total) if total else None

            dest.parent.mkdir(parents=True, exist_ok=True)
            downloaded = 0
            block_size = 8192

            with open(dest, "wb") as f:
                while True:
                    chunk = response.read(block_size)
                    if not chunk:
                        break
                    f.write(chunk)
                    downloaded += len(chunk)
                    if total:
                        pct = downloaded * 100 // total
                        mb = downloaded / (1024 * 1024)
                        total_mb = total / (1024 * 1024)
                        print(f"\r  {mb:.1f}/{total_mb:.1f} MB ({pct}%)", end="", flush=True)
                    else:
                        mb = downloaded / (1024 * 1024)
                        print(f"\r  {mb:.1f} MB downloaded", end="", flush=True)
            print()  # newline after progress
    except Exception as e:
        if dest.exists():
            dest.unlink()
        raise RuntimeError(f"Download failed: {e}") from e


def run_setup_db(args) -> int:
    """Run the database setup process."""
    project_root = Path(__file__).parent.parent.parent
    data_dir = project_root / "data"
    ictv_dir = data_dir / "ictv"
    db_path = project_root / "viroforge" / "data" / "viral_genomes.db"
    vmr_path = ictv_dir / ICTV_VMR_FILENAME

    print("=" * 60)
    print("ViroForge Database Setup")
    print("=" * 60)

    # Step 1: Download RefSeq viral genomes
    if not args.skip_download:
        print("\n[Step 1/4] Downloading RefSeq viral genomes...")
        refseq_dir = data_dir / "refseq"
        if refseq_dir.exists() and any(refseq_dir.iterdir()):
            print(f"  RefSeq data already exists at {refseq_dir}, skipping download.")
            print("  (Use --force to re-download)")
            if not args.force:
                pass  # skip
            else:
                _run_script(project_root, "scripts/download_refseq.py", ["--output", str(refseq_dir)])
        else:
            _run_script(project_root, "scripts/download_refseq.py", ["--output", str(refseq_dir)])
    else:
        print("\n[Step 1/4] Skipping RefSeq download (--download-only not needed)")

    # Step 2: Parse genomes
    if not args.download_only:
        print("\n[Step 2/4] Parsing downloaded genomes...")
        refseq_dir = data_dir / "refseq"
        parsed_dir = data_dir / "parsed"
        if not refseq_dir.exists():
            print("  ERROR: RefSeq data not found. Run without --skip-download first.")
            return 1
        _run_script(project_root, "scripts/parse_genomes.py", [
            "--input", str(refseq_dir),
            "--output", str(parsed_dir),
        ])

        # Step 3: Create and populate database
        print("\n[Step 3/4] Creating and populating database...")
        _run_script(project_root, "scripts/populate_database.py", [
            "--input", str(parsed_dir),
            "--database", str(db_path),
            "--create-db",
        ])
    else:
        print("\n[Step 2/4] Skipping parse (--download-only)")
        print("[Step 3/4] Skipping database creation (--download-only)")

    # Step 4: Download ICTV VMR and add taxonomy
    if not args.skip_taxonomy and not args.download_only:
        print("\n[Step 4/4] Adding ICTV taxonomy...")

        # Download VMR if not present
        if vmr_path.exists() and not args.force:
            print(f"  VMR file already exists at {vmr_path}, skipping download.")
        else:
            print("  Downloading ICTV Virus Metadata Resource (VMR)...")
            download_with_progress(ICTV_VMR_URL, vmr_path, label="  Fetching VMR")

        if not db_path.exists():
            print("  ERROR: Database not found. Cannot add taxonomy without database.")
            return 1

        _run_script(project_root, "scripts/parse_ictv_taxonomy.py", [
            "--vmr", str(vmr_path),
            "--output", str(ictv_dir),
            "--database", str(db_path),
            "--update-db",
        ])
    else:
        print("\n[Step 4/4] Skipping ICTV taxonomy")

    print("\n" + "=" * 60)
    if db_path.exists():
        size_mb = db_path.stat().st_size / (1024 * 1024)
        print(f"Database ready: {db_path} ({size_mb:.0f} MB)")
    print("Setup complete! You can now run: viroforge browse")
    print("=" * 60)

    return 0


def _run_script(project_root: Path, script: str, extra_args: list) -> None:
    """Run a Python script as a subprocess."""
    import subprocess
    cmd = [sys.executable, str(project_root / script)] + extra_args
    print(f"  Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, cwd=str(project_root))
    if result.returncode != 0:
        raise RuntimeError(f"Script failed with exit code {result.returncode}: {script}")
