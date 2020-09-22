.PHONY: test
test:
	poetry run pytest --cov=pycgtool

.PHONY: lint
lint:
	poetry run prospector --strictness veryhigh --tool pyflakes
# TODO: Make linting suggestions stricter once showing clean here
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
	rm -f *.itp* *.gro* *.dat* *.json* *.xtc*
	rm -rf *.ff
