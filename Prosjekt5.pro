TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

INCLUDEPATH += C:\armadillo-9.700.2\include
DEPENDPATH += C:\armadillo-9.700.2\include

LIBS += \
    -LC:\armadillo-9.700.2\examples\lib_win64 \
    -llapack_win64_MT \
    -lblas_win64_MT \

SOURCES += \
        main.cpp \
        mc_sirs.cpp \
        sirs_model.cpp

HEADERS += \
    mc_sirs.h \
    sirs_model.h
