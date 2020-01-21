# MPI Settings
QMAKE_CXX = mpicxx
QMAKE_CXX_RELEASE = $$QMAKE_CXX
QMAKE_CXX_DEBUG = $$QMAKE_CXX
QMAKE_LINK = $$QMAKE_CXX
QMAKE_CC = mpicc

QT += core
QT -= gui

CONFIG += c++11

TARGET = op_mc_hc
CONFIG += console
CONFIG -= app_bundle

TEMPLATE = app

SOURCES += main.cpp

LIBS += -L$$PWD/../partsEngine-master/ -lPartsEngine

INCLUDEPATH += $$PWD/../partsEngine-master
DEPENDPATH += $$PWD/../partsEngine-master

PRE_TARGETDEPS += $$PWD/../partsEngine-master/libPartsEngine.a

HEADERS +=
