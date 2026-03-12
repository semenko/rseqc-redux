# Contributing

## Development setup

```bash
git clone https://github.com/semenko/rseqc-redux.git
cd rseqc-redux
uv sync
```

## Running tests

```bash
# All tests
uv run pytest

# Single test file
uv run pytest tests/test_BED.py

# Single test with verbose output
uv run pytest tests/test_BED.py::test_parse_bed_line -v

# With coverage
uv run pytest --cov=rseqc --cov=scripts
```

## Linting and formatting

```bash
# Check lint
uv run ruff check .

# Auto-fix lint issues
uv run ruff check . --fix

# Format code
uv run ruff format .
```

## Type checking

```bash
uv run mypy rseqc/
```

## Building docs locally

```bash
uv sync --group docs
uv run mkdocs serve
```

Then open [http://localhost:8000](http://localhost:8000).

## Guidelines

- Write tests before changing behavior
- When fixing a bug, add a regression test first
- Target Python 3.10+
- Keep CI green: `ruff check`, `ruff format --check`, `mypy`, `pytest`
