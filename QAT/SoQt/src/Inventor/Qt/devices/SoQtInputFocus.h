#ifndef SOQT_INPUTFOCUS_H
#define SOQT_INPUTFOCUS_H

// src/Inventor/Qt/devices/SoQtInputFocus.h.  Generated from SoGuiInputFocus.h.in by configure.

/**************************************************************************\
 *
 *  This file is part of the Coin 3D visualization library.
 *  Copyright (C) by Kongsberg Oil & Gas Technologies.
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  ("GPL") version 2 as published by the Free Software Foundation.
 *  See the file LICENSE.GPL at the root directory of this source
 *  distribution for additional information about the GNU GPL.
 *
 *  For using Coin with software that can not be combined with the GNU
 *  GPL, and for taking advantage of the additional benefits of our
 *  support services, please contact Kongsberg Oil & Gas Technologies
 *  about acquiring a Coin Professional Edition License.
 *
 *  See http://www.coin3d.org/ for more information.
 *
 *  Kongsberg Oil & Gas Technologies, Bygdoy Alle 5, 0257 Oslo, NORWAY.
 *  http://www.sim.no/  sales@sim.no  coin-support@coin3d.org
 *
\**************************************************************************/

#include <Inventor/Qt/devices/SoQtDevice.h>

// *************************************************************************

class SOQT_DLL_API SoQtInputFocus : public SoQtDevice {
  SOQT_OBJECT_HEADER(SoQtInputFocus, SoQtDevice);

public:
  enum Events {
    ENTER_WINDOW = 1 << 0,
    LEAVE_WINDOW = 1 << 1,
    ALL_EVENTS   = ENTER_WINDOW | LEAVE_WINDOW
  };

  SoQtInputFocus(int mask = ALL_EVENTS);
  virtual ~SoQtInputFocus();

  virtual void enable(QWidget * widget, SoQtEventHandler * handler, void * closure);
  virtual void disable(QWidget * widget, SoQtEventHandler * handler, void * closure);

  virtual const SoEvent * translateEvent(QEvent * event);

private:
  class SoQtInputFocusP * pimpl;
  friend class SoGuiInputFocusP;
  friend class SoQtInputFocusP;
};

#define SO_QT_ALL_FOCUS_EVENTS SoQtInputFocus::ALL_EVENTS;

// *************************************************************************

#endif // ! SOQT_INPUTFOCUS_H
