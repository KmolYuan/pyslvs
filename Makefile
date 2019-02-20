# Pyslvs Makefile

# author: Yuan Chang
# copyright: Copyright (C) 2016-2019
# license: AGPL
# email: pyslvs@gmail.com

ifeq ($(OS),Windows_NT)
    SHELL = cmd
endif

all: build

.PHONY: build

# Into package folder
build:
ifeq ($(OS),Windows_NT)
	-rename __init__.py .__init__.py
	python setup.py build_ext -j0 --inplace
	-rename .__init__.py __init__.py
else
	-mv __init__.py .__init__.py
	python3 setup.py build_ext -j0 --inplace
	-mv .__init__.py __init__.py
endif

test: build
ifeq ($(OS),Windows_NT)
	python test.py
else
	python3 test.py
endif

clean:
ifeq ($(OS),Windows_NT)
	-rename .__init__.py __init__.py
	-del *.pyd /q
	-rd build /s /q
	-del src\*.c /q
	-del src\*.cpp /q
else
	-mv .__init__.py __init__.py
	-rm -f *.so
	-rm -fr build
	-rm -f src/*.c
	-rm -f src/*.cpp
endif
