TEMPLATE = app
CONFIG += console c++11
CONFIG += debug_and_release
CONFIG -= app_bundle
CONFIG -= qt

HEADERS += \
    solve.h \
    calc.h

SOURCES += \
    solve.cpp \
    derivatives.cpp \
    main.cpp \
    geometric_constraint.cpp \
    calc.cpp
