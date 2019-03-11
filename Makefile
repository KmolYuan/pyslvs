# Pyslvs Makefile

# author: Yuan Chang
# copyright: Copyright (C) 2016-2019
# license: AGPL
# email: pyslvs@gmail.com

ifeq ($(OS),Windows_NT)
    SHELL = cmd
endif

.PHONY: build test clean

all: build

# Into package folder
build: setup.py
ifeq ($(OS),Windows_NT)
	-rename __init__.py .__init__.py
	python $< build_ext -j0 --inplace
	-rename .__init__.py __init__.py
else
	-mv __init__.py .__init__.py
	python3 $< build_ext -j0 --inplace
	-mv .__init__.py __init__.py
endif

test: test_pyslvs.py build
ifeq ($(OS),Windows_NT)
	python $<
else
	python3 $<
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
