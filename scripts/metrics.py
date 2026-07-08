#!/usr/bin/env python3
"""AST-based code metrics for nh3sofc.

Counts:
  - Total function/method definitions
  - Fully annotated functions (all params + return have type hints)
  - print() calls (excluding lines with '# cli-output' comment)
  - Broad/bare except clauses (except Exception / bare except)

Usage:
    python scripts/metrics.py              # print metrics
    python scripts/metrics.py --save       # write metrics_baseline.json
    python scripts/metrics.py --check      # exit non-zero if worse than baseline
"""

from __future__ import annotations

import ast
import json
import sys
import tokenize
from pathlib import Path

PACKAGE_DIR = Path(__file__).resolve().parent.parent / "nh3sofc"
BASELINE_PATH = Path(__file__).resolve().parent.parent / "metrics_baseline.json"


def _python_files(root: Path) -> list[Path]:
    """Return all .py files under root, excluding __pycache__."""
    return sorted(
        p for p in root.rglob("*.py")
        if "__pycache__" not in str(p)
    )


def _is_fully_annotated(node: ast.FunctionDef | ast.AsyncFunctionDef) -> bool:
    """Check if a function has type annotations on all parameters and return."""
    # Check return annotation
    if node.returns is None:
        return False

    # Check all arguments (skip 'self' and 'cls')
    args = node.args
    all_args = []
    for arg_list in [args.posonlyargs, args.args, args.kwonlyargs]:
        all_args.extend(arg_list)
    if args.vararg:
        all_args.append(args.vararg)
    if args.kwarg:
        all_args.append(args.kwarg)

    for arg in all_args:
        if arg.arg in ("self", "cls"):
            continue
        if arg.annotation is None:
            return False

    return True


def _count_print_calls_tokenize(filepath: Path) -> int:
    """Count print() calls using tokenize, excluding lines with '# cli-output'."""
    count = 0
    try:
        with open(filepath, "rb") as f:
            tokens = list(tokenize.tokenize(f.readline))
    except tokenize.TokenError:
        return 0

    # Build set of line numbers that have '# cli-output' comments
    cli_output_lines: set[int] = set()
    for tok in tokens:
        if tok.type == tokenize.COMMENT and "cli-output" in tok.string:
            cli_output_lines.add(tok.start[0])

    # Find print( patterns
    for i, tok in enumerate(tokens):
        if (
            tok.type == tokenize.NAME
            and tok.string == "print"
            and tok.start[0] not in cli_output_lines
        ):
            # Check next non-NL/NEWLINE/INDENT/DEDENT/COMMENT token is '('
            for j in range(i + 1, len(tokens)):
                next_tok = tokens[j]
                if next_tok.type in (
                    tokenize.NL, tokenize.NEWLINE, tokenize.INDENT,
                    tokenize.DEDENT, tokenize.COMMENT,
                ):
                    continue
                if next_tok.type == tokenize.OP and next_tok.string == "(":
                    count += 1
                break

    return count


class _MetricsVisitor(ast.NodeVisitor):
    """AST visitor that counts functions, annotations, and broad excepts."""

    def __init__(self) -> None:
        self.total_functions: int = 0
        self.annotated_functions: int = 0
        self.broad_excepts: int = 0

    def visit_FunctionDef(self, node: ast.FunctionDef) -> None:
        self.total_functions += 1
        if _is_fully_annotated(node):
            self.annotated_functions += 1
        self.generic_visit(node)

    def visit_AsyncFunctionDef(self, node: ast.AsyncFunctionDef) -> None:
        self.total_functions += 1
        if _is_fully_annotated(node):
            self.annotated_functions += 1
        self.generic_visit(node)

    def visit_ExceptHandler(self, node: ast.ExceptHandler) -> None:
        # Bare except (no type) or except Exception
        if node.type is None:
            self.broad_excepts += 1
        elif isinstance(node.type, ast.Name) and node.type.id == "Exception":
            self.broad_excepts += 1
        self.generic_visit(node)


def collect_metrics(package_dir: Path | None = None) -> dict[str, int]:
    """Collect all metrics from the package directory."""
    if package_dir is None:
        package_dir = PACKAGE_DIR

    total_functions = 0
    annotated_functions = 0
    print_calls = 0
    broad_excepts = 0

    for filepath in _python_files(package_dir):
        # AST-based metrics
        try:
            source = filepath.read_text(encoding="utf-8")
            tree = ast.parse(source, filename=str(filepath))
        except SyntaxError:
            continue

        visitor = _MetricsVisitor()
        visitor.visit(tree)
        total_functions += visitor.total_functions
        annotated_functions += visitor.annotated_functions
        broad_excepts += visitor.broad_excepts

        # Tokenize-based print counting
        print_calls += _count_print_calls_tokenize(filepath)

    return {
        "total_functions": total_functions,
        "annotated_functions": annotated_functions,
        "print_calls": print_calls,
        "broad_excepts": broad_excepts,
    }


def save_baseline(metrics: dict[str, int], path: Path | None = None) -> None:
    """Write metrics to baseline JSON file."""
    if path is None:
        path = BASELINE_PATH
    path.write_text(json.dumps(metrics, indent=2) + "\n", encoding="utf-8")


def load_baseline(path: Path | None = None) -> dict[str, int]:
    """Load baseline metrics from JSON file."""
    if path is None:
        path = BASELINE_PATH
    return json.loads(path.read_text(encoding="utf-8"))


def check_ratchet(
    current: dict[str, int], baseline: dict[str, int]
) -> list[str]:
    """Check that metrics haven't regressed. Returns list of violations."""
    violations = []

    # These must not increase (lower is better)
    for key in ("print_calls", "broad_excepts"):
        if current[key] > baseline[key]:
            violations.append(
                f"{key}: {current[key]} > baseline {baseline[key]} (must not increase)"
            )

    # Annotation fraction must not decrease (higher is better)
    if baseline["total_functions"] > 0:
        baseline_frac = baseline["annotated_functions"] / baseline["total_functions"]
    else:
        baseline_frac = 0.0

    if current["total_functions"] > 0:
        current_frac = current["annotated_functions"] / current["total_functions"]
    else:
        current_frac = 0.0

    if current_frac < baseline_frac - 1e-9:  # small tolerance for rounding
        violations.append(
            f"annotation_fraction: {current_frac:.4f} < baseline {baseline_frac:.4f} "
            f"(must not decrease)"
        )

    return violations


def main() -> int:
    metrics = collect_metrics()

    print(f"total_functions:     {metrics['total_functions']}")
    print(f"annotated_functions: {metrics['annotated_functions']}")
    annotation_pct = (
        100 * metrics["annotated_functions"] / metrics["total_functions"]
        if metrics["total_functions"] > 0
        else 0
    )
    print(f"annotation_fraction: {annotation_pct:.1f}%")
    print(f"print_calls:         {metrics['print_calls']}")
    print(f"broad_excepts:       {metrics['broad_excepts']}")

    if "--save" in sys.argv:
        save_baseline(metrics)
        print(f"\nBaseline saved to {BASELINE_PATH}")

    if "--check" in sys.argv:
        if not BASELINE_PATH.exists():
            print(f"\nNo baseline found at {BASELINE_PATH}; run with --save first.")
            return 1
        baseline = load_baseline()
        violations = check_ratchet(metrics, baseline)
        if violations:
            print("\nRatchet violations:")
            for v in violations:
                print(f"  FAIL: {v}")
            return 1
        else:
            print("\nAll ratchet checks passed.")

    return 0


if __name__ == "__main__":
    sys.exit(main())
