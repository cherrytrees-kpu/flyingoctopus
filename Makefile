# Run program
run:
	python VariantAnnotation.py

# Install dependencies
install:
	pip install -r requirements.txt

# Freeze current dependencies to requirements file
freeze:
	pip freeze > requirements.txt

.PHONY: run install freeze

# Virtual environment specific commands

# Ensure user has virtualenv
venv-requirements:
	virtualenv --version || (echo "please install virtualenv" && exit 1)

# Generate new virtual environment
venv-fresh:
	virtualenv venv

# Create a fresh venv for this project with dependencies installed
venv-install: venv-requirements venv-fresh
	( \
		source venv/bin/activate; \
		pip install -r requirements.txt; \
	)

# Freeze venv (can also run `make freeze` while venv is active instead)
venv-freeze: venv-requirements
	( \
		source venv/bin/activate; \
		pip freeze > requirements.txt; \
	)

.PHONY: venv-install venv-activate venv-freeze