"""Ratchet test: code metrics must not regress below baseline.

This test enforces monotonic improvement:
  - print_calls and broad_excepts must not increase
  - annotation fraction must not decrease

To update the baseline after cleaning a module, run:
    python scripts/metrics.py --save

Then commit the updated metrics_baseline.json alongside your cleanup.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

# Import from the metrics script
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))
from metrics import collect_metrics, check_ratchet  # noqa: E402

BASELINE_PATH = Path(__file__).resolve().parent.parent / "metrics_baseline.json"


@pytest.fixture
def baseline() -> dict[str, int]:
    if not BASELINE_PATH.exists():
        pytest.skip("No metrics_baseline.json found; run 'python scripts/metrics.py --save'")
    return json.loads(BASELINE_PATH.read_text(encoding="utf-8"))


@pytest.fixture
def current_metrics() -> dict[str, int]:
    return collect_metrics()


class TestMetricsRatchet:
    """Ensure code quality metrics never regress."""

    def test_print_calls_not_increased(
        self, current_metrics: dict[str, int], baseline: dict[str, int]
    ) -> None:
        assert current_metrics["print_calls"] <= baseline["print_calls"], (
            f"print() calls increased: {current_metrics['print_calls']} > "
            f"baseline {baseline['print_calls']}. "
            f"Use logging instead, or add '# cli-output' comment if intentional."
        )

    def test_broad_excepts_not_increased(
        self, current_metrics: dict[str, int], baseline: dict[str, int]
    ) -> None:
        assert current_metrics["broad_excepts"] <= baseline["broad_excepts"], (
            f"Broad excepts increased: {current_metrics['broad_excepts']} > "
            f"baseline {baseline['broad_excepts']}. "
            f"Use specific exception types."
        )

    def test_annotation_fraction_not_decreased(
        self, current_metrics: dict[str, int], baseline: dict[str, int]
    ) -> None:
        if baseline["total_functions"] == 0:
            return

        baseline_frac = baseline["annotated_functions"] / baseline["total_functions"]
        current_frac = (
            current_metrics["annotated_functions"] / current_metrics["total_functions"]
            if current_metrics["total_functions"] > 0
            else 0.0
        )
        assert current_frac >= baseline_frac - 1e-9, (
            f"Annotation fraction decreased: {current_frac:.4f} < "
            f"baseline {baseline_frac:.4f}. "
            f"Add type annotations to new functions."
        )

    def test_no_ratchet_violations(
        self, current_metrics: dict[str, int], baseline: dict[str, int]
    ) -> None:
        """Comprehensive check using the same logic as scripts/metrics.py --check."""
        violations = check_ratchet(current_metrics, baseline)
        assert not violations, (
            "Ratchet violations:\n" + "\n".join(f"  - {v}" for v in violations)
        )
