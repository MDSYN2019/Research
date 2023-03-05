VERSION=2.0.0

isEmpty(PREFIX) {
  PREFIX=usr/local
}
  
isEmpty(QATLIBDIR) {
  QATLIBDIR=/$$PREFIX/lib
}

TEMPLATE = lib dll
TARGET = QatDataAnalysis
DEPENDPATH += .
INCLUDEPATH += . .. ../../QatGenericFunctions

DESTDIR=../../../lib
CONFIG += build_all release c++11


QMAKE_CXXFLAGS += -DQATLIBDIR=$$QATLIBDIR
# Input
HEADERS +=            ../QatDataAnalysis/*.h ../QatDataAnalysis/*.icc

SOURCES +=            *.cpp

target.path=/$$QATLIBDIR
INSTALLS += target

headers.path=/$$PREFIX/include/QatDataAnalysis
headers.files= ../QatDataAnalysis/*h ../QatDataAnalysis/*.icc
INSTALLS += headers

pc.path  = $$QATLIBDIR/pkgconfig
mac {
 pc.files = ../pkgconfig/mac/*.pc
}
linux {
 pc.files = ../pkgconfig/linux/*.pc
}
INSTALLS += pc

LIBS += -L$(DESTDIR) -lQatGenericFunctions -ldl


