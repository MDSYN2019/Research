#ifndef COIN_QTCONFIG_H
#define COIN_QTCONFIG_H

/* Define this if QApplication::hasPendingEvents() is available */
#define HAVE_QAPPLICATION_HASPENDINGEVENTS 1

/* Define this to 1 if operator==(QGLFormat&, QGLFormat&) is available */
#define HAVE_QGLFORMAT_EQ_OP 1

/* Define this to 1 if QGLFormat::setOverlay() is available */
#define HAVE_QGLFORMAT_SETOVERLAY 1

/* Define this to 1 if QGLWidget::setAutoBufferSwap() is available */
#define HAVE_QGLWIDGET_SETAUTOBUFFERSWAP 1

/* Define to 1 if you have the <qstylefactory.h> header file. */
/** QStyleFactory was added in Qt 3.0. **/
#define HAVE_QSTYLEFACTORY_H 1

/* Define this if QWidget::showFullScreen() is available */
#define HAVE_QWIDGET_SHOWFULLSCREEN 1

#define HAVE_QWIDGET_SETWINDOWSTATE 1



/* QGLFormat::setSampleBuffers() was introduced in Qt 4.0 */
#define HAVE_QGLFORMAT_SETSAMPLEBUFFERS 1

#endif
