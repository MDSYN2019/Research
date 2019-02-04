VERSION = 2.0.0

isEmpty(PREFIX) {
  PREFIX=/usr/local
}
  
isEmpty(QATLIBDIR) {
  QATLIBDIR=$$PREFIX/lib
}


TEMPLATE = lib dll
TARGET = SoQt
DEPENDPATH += .  devices editors engines nodes viewers widgets 
INCLUDEPATH += . devices editors
INCLUDEPATH += . ../.. /usr/local/include
LIBS +=  -L$$QATLIBDIR -lCoin

CONFIG += qt release c++11
QT     += widgets opengl


QMAKE_CXXFLAGS += -DSOQT_INTERNAL -DSOQT_DEBUG=1 -DQT_SHARED -DNDEBUG

QMAKE_CXXFLAGS +=  -fno-builtin 

# Input
HEADERS += SoAny.h \
           SoGuiComponentP.h \
           SoGuiGLWidgetP.h \
           SoGuiP.h \
           SoQt.h \
           SoQtBasic.h \
           SoQtColorEditor.h \
           SoQtComponent.h \
           SoQtComponentP.h \
           SoQtCursor.h \
           SoQtGLWidget.h \
           SoQtGLWidgetP.h \
           SoQtImageReader.h \
           SoQtInternal.h \
           SoQtLightSliderSet.h \
           SoQtMaterialEditor.h \
           SoQtMaterialSliderSet.h \
           SoQtObject.h \
           SoQtP.h \
           SoQtRenderArea.h \
           SoQtSignalThread.h \
           SoQtSliderSet.h \
           SoQtSliderSetBase.h \
           SoQtTransformSliderSet.h \
           devices/6DOFEvents.h \
           devices/SoGuiDeviceP.h \
           devices/SoGuiInputFocusP.h \
           devices/SoGuiKeyboardP.h \
           devices/SoGuiMouseP.h \
           devices/SoQtDevice.h \
           devices/SoQtDeviceP.h \
           devices/SoQtInputFocus.h \
           devices/SoQtKeyboard.h \
           devices/SoQtMouse.h \
           devices/spwinput.h \
           devices/spwinput_win32.h \
           editors/RadioGroupKit.h \
           editors/RGBCubeEditorKit.h \
           editors/SoQtColorEditor.h \
           editors/SoQtMaterialEditor.h \
           engines/SoGuiEngines.h \
           engines/SoGuiFormat.h \
           engines/SoGuiRadioGroup.h \
           nodes/SoGuiClickCounter.h \
           nodes/SoGuiColorEditor.h \
           nodes/SoGuiFrame.h \
           nodes/SoGuiImage.h \
           nodes/SoGuiLabel.h \
           nodes/SoGuiMaterialEditor.h \
           nodes/SoGuiNodes.h \
           nodes/SoGuiPane.h \
           nodes/SoGuiPosition.h \
           nodes/SoGuiRadioButton.h \
           nodes/SoGuiSceneTexture2.h \
           nodes/SoGuiSlider1.h \
           nodes/SoGuiSlider2.h \
           nodes/SoGuiToggleButton.h \
           nodes/SoGuiTranslation.h \
           nodes/SoGuiViewpointWrapper.h \
           nodes/SoGuiViewportFix.h \
           viewers/SoGuiExaminerViewerP.h \
           viewers/SoGuiFullViewerP.h \
           viewers/SoGuiPlaneViewerP.h \
           viewers/SoGuiViewerP.h \
           viewers/SoQtConstrainedViewer.h \
           viewers/SoQtExaminerViewer.h \
           viewers/SoQtExaminerViewerP.h \
           viewers/SoQtFlyViewer.h \
           viewers/SoQtFullViewer.h \
           viewers/SoQtFullViewerP.h \
           viewers/SoQtPlaneViewer.h \
           viewers/SoQtPlaneViewerP.h \
           viewers/SoQtViewer.h \
           widgets/QtNativePopupMenu.h \
           widgets/SoAnyThumbWheel.h \
           widgets/SoQtGLArea.h \
           widgets/SoQtPopupMenu.h \
           widgets/SoQtThumbWheel.h 
SOURCES += SoAny.cpp \
           SoQt.cpp \
           SoQtCommon.cpp \
           SoQtComponent.cpp \
           SoQtComponentCommon.cpp \
           SoQtCursor.cpp \
           SoQtGLWidget.cpp \
           SoQtGLWidgetCommon.cpp \
           SoQtImageReader.cpp \
           SoQtLightSliderSet.cpp \
           SoQtObject.cpp \
           SoQtRenderArea.cpp \
           SoQtSignalThread.cpp \
           SoQtSliderSet.cpp \
           SoQtSliderSetBase.cpp \
           SoQtTransformSliderSet.cpp \
           devices/6DOFEvents.cpp \
           devices/SoQtDevice.cpp \
           devices/SoQtDeviceCommon.cpp \
           devices/SoQtInputFocus.cpp \
           devices/SoQtInputFocusCommon.cpp \
           devices/SoQtKeyboard.cpp \
           devices/SoQtKeyboardCommon.cpp \
           devices/SoQtMouse.cpp \
           devices/SoQtMouseCommon.cpp \
           devices/spwinput_win32.c \
           devices/spwinput_x11.cpp \
           editors/RadioGroupKit.cpp \
           editors/RGBCubeEditorKit.cpp \
           editors/SoQtColorEditor.cpp \
           editors/SoQtMaterialEditor.cpp \
           engines/Engines.cpp \
           engines/Format.cpp \
           engines/RadioGroup.cpp \
           nodes/ClickCounter.cpp \
           nodes/ColorEditor.cpp \
           nodes/Frame.cpp \
           nodes/Image.cpp \
           nodes/Label.cpp \
           nodes/MaterialEditor.cpp \
           nodes/Nodes.cpp \
           nodes/Pane.cpp \
           nodes/Position.cpp \
           nodes/RadioButton.cpp \
           nodes/SceneTexture2.cpp \
           nodes/Slider1.cpp \
           nodes/Slider2.cpp \
           nodes/ToggleButton.cpp \
           nodes/Translation.cpp \
           nodes/ViewpointWrapper.cpp \
           nodes/ViewportFix.cpp \
           viewers/ExaminerViewer.cpp \
           viewers/FullViewer.cpp \
           viewers/PlaneViewer.cpp \
           viewers/SoQtConstrainedViewer.cpp \
           viewers/SoQtExaminerViewer.cpp \
           viewers/SoQtFlyViewer.cpp \
           viewers/SoQtFullViewer.cpp \
           viewers/SoQtPlaneViewer.cpp \
           viewers/SoQtViewer.cpp \
           widgets/QtNativePopupMenu.cpp \
           widgets/SoAnyThumbWheel.cpp \
           widgets/SoQtGLArea.cpp \
           widgets/SoQtPopupMenu.cpp \
           widgets/SoQtThumbWheel.cpp 


target.path=$$QATLIBDIR
INSTALLS += target



INSTALL_PREFIX = $$PREFIX/include/Inventor/Qt
INSTALL_HEADERS =  \
           SoAny.h \
           SoQt.h \
           SoQtBasic.h \
           SoQtColorEditor.h \
           SoQtComponent.h \
           SoQtCursor.h \
           SoQtGLWidget.h \
           SoQtImageReader.h \
           SoQtInternal.h \
           SoQtLightSliderSet.h \
           SoQtMaterialEditor.h \
           SoQtMaterialSliderSet.h \
           SoQtObject.h \
           SoQtRenderArea.h \
           SoQtSignalThread.h \
           SoQtSliderSet.h \
           SoQtSliderSetBase.h \
           SoQtTransformSliderSet.h \
           devices/6DOFEvents.h \
           devices/SoQtDevice.h \
           devices/SoQtInputFocus.h \
           devices/SoQtKeyboard.h \
           devices/SoQtMouse.h \
           devices/spwinput.h \
           devices/spwinput_win32.h \
           editors/RadioGroupKit.h \
           editors/RGBCubeEditorKit.h \
           editors/SoQtColorEditor.h \
           editors/SoQtMaterialEditor.h \
           engines/SoGuiEngines.h \
           engines/SoGuiFormat.h \
           engines/SoGuiRadioGroup.h \
           nodes/SoGuiClickCounter.h \
           nodes/SoGuiColorEditor.h \
           nodes/SoGuiFrame.h \
           nodes/SoGuiImage.h \
           nodes/SoGuiLabel.h \
           nodes/SoGuiMaterialEditor.h \
           nodes/SoGuiNodes.h \
           nodes/SoGuiPane.h \
           nodes/SoGuiPosition.h \
           nodes/SoGuiRadioButton.h \
           nodes/SoGuiSceneTexture2.h \
           nodes/SoGuiSlider1.h \
           nodes/SoGuiSlider2.h \
           nodes/SoGuiToggleButton.h \
           nodes/SoGuiTranslation.h \
           nodes/SoGuiViewpointWrapper.h \
           nodes/SoGuiViewportFix.h \
           viewers/SoQtConstrainedViewer.h \
           viewers/SoQtExaminerViewer.h \
           viewers/SoQtFlyViewer.h \
           viewers/SoQtFullViewer.h \
           viewers/SoQtPlaneViewer.h \
           viewers/SoQtViewer.h \
           widgets/QtNativePopupMenu.h \
           widgets/SoAnyThumbWheel.h \
           widgets/SoQtGLArea.h \
           widgets/SoQtPopupMenu.h \
           widgets/SoQtThumbWheel.h 

           
include(headerinstall.pri)

pc.path  = $$QATLIBDIR/pkgconfig
mac {
    pc.files = ../../../pkgconfig/mac/*.pc
}
linux {
    pc.files = ../../../pkgconfig/linux/*.pc
}  
INSTALLS += pc

mac {
  CONFIG -= app_bundle
}

