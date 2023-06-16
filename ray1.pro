QT += gui widgets concurrent

CONFIG += c++1z optimize_full
CONFIG -= app_bundle

# You can make your code fail to compile if it uses deprecated APIs.
# In order to do so, uncomment the following line.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

SOURCES += \
	bvh_node.cpp \
        main.cpp \
	rayview_rtiow.cpp \
        rayview.cpp \
        hitPosition.cpp \
	calc.cpp

# Default rules for deployment.
qnx: target.path = /tmp/$${TARGET}/bin
else: unix:!android: target.path = /opt/$${TARGET}/bin
!isEmpty(target.path): INSTALLS += target

FORMS += \
	rayview_rtiow.ui \
	rayview.ui

HEADERS += \
	Ray.h \
	bvh_node.h \
	calc.h \
	camera.h \
	hitPosition.h \
	rayview_rtiow.h \
	rayview.h \
	shapes.h
	
