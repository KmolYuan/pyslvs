TEMPLATE = app
CONFIG += console c++11
CONFIG += debug_and_release
CONFIG -= app_bundle
CONFIG -= qt

HEADERS += solve.h \
    constrain_func.h \
    position.h

SOURCES += \
    solve.cpp \
    constraint_func.cpp \
    derivatives.cpp \
    main.cpp \
    position.cpp
