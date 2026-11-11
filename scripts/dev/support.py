#!/usr/bin/env python3
"""Print the SPEC 0 support windows for this project."""

from __future__ import annotations

import argparse
import json
import re
import shlex
import subprocess
import tomllib
from datetime import date, datetime, timedelta
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

from packaging.requirements import Requirement
from packaging.specifiers import InvalidSpecifier, SpecifierSet
from packaging.version import InvalidVersion, Version
from rich.console import Console
from rich.table import Table

PYTHON_RELEASES_URL = "https://endoflife.date/api/python.json"
PYTHON_SUPPORT = timedelta(days=365 * 3)
DEPENDENCY_SUPPORT = timedelta(days=365 * 2)
PROJECT_FILE = Path(__file__).resolve().parents[2] / "pyproject.toml"
console = Console()


def satisfies_python(version: str, requirement: str | None) -> bool:
    """Evaluate a PyPI Requires-Python value."""

    if not requirement:
        return True

    try:
        return SpecifierSet(requirement).contains(Version(version))
    except (InvalidSpecifier, InvalidVersion):
        return False


def project_dependencies() -> list[tuple[str, str | None]]:
    with PROJECT_FILE.open("rb") as stream:
        project = tomllib.load(stream)["project"]

    dependencies = [(requirement, None) for requirement in project.get("dependencies", [])]
    for group, requirements in project.get("optional-dependencies", {}).items():
        dependencies.extend((requirement, group) for requirement in requirements)

    packages = {(Requirement(requirement).name, group) for requirement, group in dependencies}
    return sorted(packages, key=lambda package: package[0].casefold())


def fetch_json(url: str, label: str) -> dict | list:
    request = Request(url, headers={"Accept": "application/json"})
    try:
        with urlopen(request, timeout=30) as response:
            return json.load(response)
    except (HTTPError, URLError) as error:
        raise RuntimeError(f"could not query {label}: {error}") from error


def python_releases() -> dict[str, date]:
    releases = {}
    for release in fetch_json(PYTHON_RELEASES_URL, "the Python release calendar"):
        cycle = release.get("cycle")
        release_date = release.get("releaseDate")
        if not isinstance(cycle, str) or not re.fullmatch(r"\d+\.\d+", cycle) or not release_date:
            continue
        releases[cycle] = date.fromisoformat(release_date)
    return releases


def fetch_package(name: str) -> dict:
    return fetch_json(f"https://pypi.org/pypi/{name}/json", f"PyPI package {name}")


def release_date(files: list[dict]) -> date | None:
    timestamps = [file["upload_time_iso_8601"] for file in files if file.get("upload_time_iso_8601")]
    if not timestamps:
        return None
    return min(datetime.fromisoformat(timestamp.replace("Z", "+00:00")) for timestamp in timestamps).date()


def supported_releases(name: str, python_version: str, today: date) -> list[tuple[str, date]]:
    package = fetch_package(name)
    candidates = []
    for version, files in package.get("releases", {}).items():
        try:
            parsed_version = Version(version)
        except InvalidVersion:
            continue
        feature_version = parsed_version.release + (0,) * (3 - len(parsed_version.release))
        if parsed_version.is_prerelease or parsed_version.is_devrelease or feature_version[2] != 0:
            continue
        usable_files = [file for file in files if not file.get("yanked")]
        if not usable_files or not any(
            satisfies_python(python_version, file.get("requires_python")) for file in usable_files
        ):
            continue
        published = release_date(usable_files)
        if published is not None and published + DEPENDENCY_SUPPORT >= today:
            candidates.append((version, published))

    return sorted(candidates, key=lambda candidate: Version(candidate[0]))


def update_dependencies(updates: dict[str | None, list[str]], dry_run: bool = False) -> None:
    commands = []
    for group, requirements in updates.items():
        command = ["uv", "add"]
        if group is not None:
            command.extend(["--optional", group])
        command.extend(requirements)
        commands.append(command)

    for command in commands:
        label = "Would run:" if dry_run else "Running:"
        console.print(f"[dim]{label}[/dim] {shlex.join(command)}")
        if not dry_run:
            subprocess.run(command, check=True)


def update_python_requirement(version: str, dry_run: bool = False) -> None:
    contents = PROJECT_FILE.read_text()
    updated, replacements = re.subn(
        r'^requires-python\s*=\s*["\'][^"\']*["\']$',
        f'requires-python = ">={version}"',
        contents,
        count=1,
        flags=re.MULTILINE,
    )
    if replacements != 1:
        raise RuntimeError(f"could not find requires-python in {PROJECT_FILE}")

    label = "Would update" if dry_run else "Updating"
    console.print(f"[dim]{label}:[/dim] {PROJECT_FILE} requires-python >= {version}")
    if not dry_run:
        PROJECT_FILE.write_text(updated)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--update",
        action="store_true",
        help="update the Python minimum and dependency minimums",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="show planned metadata and uv add changes without applying them",
    )
    args = parser.parse_args()

    today = date.today()
    supported_python = [
        (version, released, released + PYTHON_SUPPORT)
        for version, released in sorted(python_releases().items(), key=lambda item: Version(item[0]))
        if released + PYTHON_SUPPORT >= today
    ]
    if not supported_python:
        console.print("No Python release is currently within the SPEC 0 support window.", style="red")
        return 1

    oldest_version = min(supported_python, key=lambda item: Version(item[0]))[0]
    python_table = Table(title="Python Support Windows (SPEC 0: 3 years)")
    python_table.add_column("Version", style="cyan")
    python_table.add_column("Released")
    python_table.add_column("Supported until")
    for version, released, drop_date in supported_python:
        marker = " (oldest)" if version == oldest_version else ""
        python_table.add_row(f"{version}{marker}", str(released), str(drop_date))
    console.print(python_table)

    updates: dict[str | None, list[str]] = {}
    dependency_rows = []
    for name, group in project_dependencies():
        releases = supported_releases(name, oldest_version, today)
        if not releases:
            dependency_rows.append((name, "none", "none", "no compatible release"))
            continue
        oldest, newest = releases[0], releases[-1]
        updates.setdefault(group, []).append(f"{name}>={oldest[0]}")
        newest_version = newest[0] if newest != oldest else "-"
        window = f"{oldest[1]} -> {oldest[1] + DEPENDENCY_SUPPORT}"
        dependency_rows.append((name, oldest[0], newest_version, window))

    dependency_table = Table(title="Dependency Support Windows (SPEC 0: 2 years)")
    dependency_table.add_column("Package", style="cyan")
    dependency_table.add_column("Oldest supported", style="green")
    dependency_table.add_column("Newest supported", style="green")
    dependency_table.add_column("Oldest release window")
    for row in dependency_rows:
        dependency_table.add_row(*row)
    console.print(dependency_table)

    if args.update or args.dry_run:
        update_python_requirement(oldest_version, dry_run=args.dry_run)
        update_dependencies(updates, dry_run=args.dry_run)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
