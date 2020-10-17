.PHONY: test
test:
# Collect coverage data, but don't report it - for CI build
	poetry run pytest --cov=pycgtool --cov-report=

.PHONY: cov
cov:
	poetry run pytest --cov=pycgtool --cov-fail-under=80

.PHONY: lint
lint:
# Lint using only pyflakes - checks for actual errors
	poetry run prospector --strictness veryhigh --test-warnings --tool pyflakes
# Lint using range of tools, but don't fail the build if we get warnings
	poetry run prospector --strictness veryhigh --test-warnings --member-warnings --max-line-length 120 --zero-exit

.PHONY: docs
docs:
	SPHINXBUILD="poetry run sphinx-build" make -C docs html

.PHONY: build
build:
	poetry build

.PHONY: publish
publish:
	rm -rf dist/
	poetry version prerelease
	poetry build
	poetry publish -r testpypi

.PHONY: clean
clean:
	rm -f .coverage *.itp* *.gro* *.dat* *.json* *.xtc*
	rm -rf *.ff
