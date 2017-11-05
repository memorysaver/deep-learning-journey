default: sync-local

sync-local: compile-deps
		pip-sync requirements.txt

compile-deps:
		pip-compile requirements/local.ini

install: sync-local
		pip install -r requirements.txt
		pip install matplotlib
