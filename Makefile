.PHONY: test
test:
	poetry run pytest --cov=pycgtool

.PHONY: lint
lint:
	poetry run prospector --tool pyflakes
# TODO: Make linting suggestions stricter once showing clean here
	poetry run prospector --strictness medium --test-warnings --member-warnings --max-line-length 120 --zero-exit

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
