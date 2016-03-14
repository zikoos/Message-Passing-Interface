TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.c

HEADERS += \
    rasterfile.h

LIBS += -L/usr/lib/openmpi/
INCLUDEPATH += /usr/include/openmpi/ompi/mpi/cxx

