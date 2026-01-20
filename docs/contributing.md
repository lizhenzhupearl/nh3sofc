# Contributing

Thank you for your interest in contributing to NH3-SOFC!

## Development Setup

1. Clone the repository:

```bash
git clone https://github.com/yourusername/NH3SOFC.git
cd NH3SOFC
```

2. Create a development environment:

```bash
conda create -n nh3sofc-dev python=3.10
conda activate nh3sofc-dev
pip install -e ".[dev]"
```

3. Install pre-commit hooks:

```bash
pre-commit install
```

## Code Style

- We use [Black](https://black.readthedocs.io/) for code formatting
- We use [isort](https://pycqa.github.io/isort/) for import sorting
- We use [flake8](https://flake8.pycqa.org/) for linting

Run formatting:

```bash
black nh3sofc/
isort nh3sofc/
flake8 nh3sofc/
```

## Testing

Run tests with pytest:

```bash
pytest tests/
pytest tests/ -v --cov=nh3sofc
```

## Pull Request Process

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Make your changes
4. Run tests and ensure they pass
5. Commit your changes (`git commit -m 'Add amazing feature'`)
6. Push to the branch (`git push origin feature/amazing-feature`)
7. Open a Pull Request

## Reporting Issues

Please use the [GitHub Issues](https://github.com/yourusername/NH3SOFC/issues) page to report bugs or request features.

When reporting bugs, please include:

- Python version
- ASE version
- Operating system
- Complete error traceback
- Minimal code to reproduce the issue

## Documentation

Documentation is built with MkDocs:

```bash
pip install mkdocs mkdocs-material mkdocstrings[python]
mkdocs serve  # Local preview at http://127.0.0.1:8000
mkdocs build  # Build static site
```
