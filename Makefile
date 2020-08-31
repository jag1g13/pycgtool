.PHONY: test
test:
	poetry run pytest --cov=pycgtool

.PHONY: lint
lint:
	poetry run prospector --tool pyflakes
	poetry run prospector --strictness veryhigh --zero-exit

.PHONY: build
build:
	poetry build

.PHONY: publish
publish:
	rm -rf dist/
	poetry version prerelease
	poetry build
	poetry publish -r testpypi
