import nox

nox.options.default_venv_backend = "uv"


@nox.session(python=["3.10", "3.11", "3.12", "3.13"])
def tests(session: nox.Session) -> None:
    """Run the test suite."""
    session.install(".[dev]")
    session.run("pytest", "--cov=rseqc", "--cov-report=xml", *session.posargs)


@nox.session
def lint(session: nox.Session) -> None:
    """Run ruff linter."""
    session.install("ruff")
    session.run("ruff", "check", ".")


@nox.session
def format_check(session: nox.Session) -> None:
    """Check code formatting with ruff."""
    session.install("ruff")
    session.run("ruff", "format", "--check", ".")


@nox.session
def typecheck(session: nox.Session) -> None:
    """Run mypy type checker."""
    session.install(".[dev]")
    session.run("mypy", "rseqc/")
