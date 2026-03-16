"""Command-line entrypoint for the mutation generation pipeline."""

from pathlib import Path
import sys

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

from src.pipeline import run_pipeline


if __name__ == "__main__":
    run_pipeline()
    print("Pipeline completed. Mutant structures written to results/structures")
