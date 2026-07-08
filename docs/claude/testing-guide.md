# Testing Guide

## Run Tests Before and After Changes

```bash
pytest tests/            # All tests
pytest tests/ -v         # Verbose output
pytest tests/test_regressions.py  # Specific file
```

## Test Structure

| File | Purpose |
|------|---------|
| `tests/conftest.py` | Shared fixtures (mock structures) |
| `tests/test_smoke.py` | Basic imports and instantiation |
| `tests/test_regressions.py` | Tests for bugs that were fixed |

## When Developing

1. **After making changes:** Run `pytest tests/` to catch regressions
2. **When fixing a bug:** Add a regression test to `tests/test_regressions.py`
3. **When adding features:** Add smoke tests for new classes/methods

## Adding Regression Tests

When fixing a bug, add a test to `tests/test_regressions.py`:

```python
class TestYourBugFix:
    """
    Bug: [Description of the bug]
    Fix: [How it was fixed]
    Date fixed: YYYY-MM-DD
    """

    def test_bug_is_fixed(self, fixture):
        # Test that verifies the fix works
        assert result == expected
```

## CI/CD

GitHub Actions runs tests automatically on push/PR to `main`. See `.github/workflows/ci.yml`.
