// src/Inventor/Qt/devices/SoQtKeyboardCommon.cpp.  Generated from SoGuiKeyboard.cpp.in by configure.

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

/*!
  \class SoQtKeyboard SoQtKeyboard.h Inventor/Qt/devices/SoQtKeyboard.h
  \brief The SoQtKeyboard class is the keyboard input device abstraction.
  \ingroup devices

  The SoQtKeyboard class is the glue between native keyboard
  handling and keyboard interaction with the Inventor scenegraph.

  All components derived from the SoQtRenderArea have got an
  SoQtKeyboard device attached by default.
*/

// *************************************************************************

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif // HAVE_CONFIG_H

#include <Inventor/Qt/devices/SoQtKeyboard.h>
#include <Inventor/Qt/devices/SoGuiKeyboardP.h>
#include <Inventor/events/SoKeyboardEvent.h>

// *************************************************************************

SOQT_OBJECT_SOURCE(SoQtKeyboard);

// *************************************************************************

/*!
  \enum SoQtKeyboard::Events
  Enumeration over supported event types.
*/

/*!
  \var SoQtKeyboard::Events SoQtKeyboard::KEY_PRESS
  Maskbit for a keyboard button press event.
*/

/*!
  \var SoQtKeyboard::Events SoQtKeyboard::KEY_RELEASE
  Maskbit for a keyboard button release event.
*/

/*!
  \var SoQtKeyboard::Events SoQtKeyboard::ALL_EVENTS
  Combined bitmask for all possible events.
*/

/*!
  \fn SoQtKeyboard::SoQtKeyboard(int mask)

  Constructor. The \a mask specifies which keyboard-related events to
  handle. Others will just be ignored.
*/

/*!
  \fn SoQtKeyboard::~SoQtKeyboard()

  Destructor.
*/

// *************************************************************************

#ifndef DOXYGEN_SKIP_THIS

SoGuiKeyboardP::SoGuiKeyboardP(void)
{
  this->kbdevent = new SoKeyboardEvent;
}

SoGuiKeyboardP::~SoGuiKeyboardP()
{
  delete this->kbdevent;
}

#endif // !DOXYGEN_SKIP_THIS

// *************************************************************************
