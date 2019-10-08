# Set Configuration Variables
NAME = pdbt

# Configure Verbose Console Output
# If V is defined (by export $V=True or adding V=TRUE to the make command), make will display verbose command information.
ifdef V
	Q =
else
# Commands proceeded with an @ will not be returned to the console.
	Q = @
endif

#### End of System Configuration ####

.PHONY: init build clean distclean install install-dev test uninstall

init:
	@echo "Initializing..."
	@echo "Installing Required Packages..."
	pip install -r requirements.txt
	@echo "Done."

build:
	@echo "Building Source Distribution..."
	python setup.py sdist
	@echo "Building Binary Distribution..."
	python setup.py bdist_egg
	python setup.py bdist_wheel
	@echo "Done."

clean:

distclean:
	@echo "Cleaning..."
	@echo "Removing build and dist folders..."
	$(Q)rm -rf dist build
	@echo "Remove setuptools folders..."
	$(Q)rm -rf pdbt.egg-info/
	$(Q)rm -rf pdbt/pdbt.egg-info/
	@echo "Done."

install:
	@echo "Installing PDBT from Github Repository..."
	$(Q)pip install git+https://github.com/gillejw/PDBT.git#egg=pdbt

install-dev:
	@echo "Installing PDBT from Local Source..."
	$(Q)pip install -e .[dev]
	@echo "Done."

lint:
	python setup.py flake8

test:
	@echo "Running Tests..."
	pytest --cache-clear
	@echo "Done."

test2:
	python setup.py test

uninstall:
	@echo "Uninstalling PDBT..."
	pip uninstall pdbt -y
	@echo "Done."
