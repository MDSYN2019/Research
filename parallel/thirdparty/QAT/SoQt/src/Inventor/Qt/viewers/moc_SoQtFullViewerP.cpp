/****************************************************************************
** Meta object code from reading C++ file 'SoQtFullViewerP.h'
**
** Created by: The Qt Meta Object Compiler version 63 (Qt 4.8.6)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "SoQtFullViewerP.h"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'SoQtFullViewerP.h' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 63
#error "This file was generated using the moc from 4.8.6. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_SoQtFullViewerP[] = {

 // content:
       6,       // revision
       0,       // classname
       0,    0, // classinfo
      20,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       0,       // signalCount

 // slots: signature, parameters, type, tag, flags
      17,   16,   16,   16, 0x0a,
      42,   36,   16,   16, 0x0a,
      66,   16,   16,   16, 0x0a,
      86,   16,   16,   16, 0x0a,
     106,   36,   16,   16, 0x0a,
     131,   16,   16,   16, 0x0a,
     152,   16,   16,   16, 0x0a,
     173,   36,   16,   16, 0x0a,
     199,   16,   16,   16, 0x0a,
     221,   16,   16,   16, 0x0a,
     245,   16,   16,   16, 0x0a,
     265,   16,   16,   16, 0x0a,
     285,   16,   16,   16, 0x0a,
     308,   16,   16,   16, 0x0a,
     331,   16,   16,   16, 0x0a,
     351,   16,   16,   16, 0x0a,
     369,   16,   16,   16, 0x0a,
     390,   16,   16,   16, 0x0a,
     410,   16,   16,   16, 0x0a,
     437,   16,   16,   16, 0x0a,

       0        // eod
};

static const char qt_meta_stringdata_SoQtFullViewerP[] = {
    "SoQtFullViewerP\0\0leftWheelPressed()\0"
    "value\0leftWheelChanged(float)\0"
    "leftWheelReleased()\0rightWheelPressed()\0"
    "rightWheelChanged(float)\0rightWheelReleased()\0"
    "bottomWheelPressed()\0bottomWheelChanged(float)\0"
    "bottomWheelReleased()\0interactbuttonClicked()\0"
    "viewbuttonClicked()\0homebuttonClicked()\0"
    "sethomebuttonClicked()\0viewallbuttonClicked()\0"
    "seekbuttonClicked()\0selectedViewing()\0"
    "selectedDecoration()\0selectedHeadlight()\0"
    "increaseInteractiveCount()\0"
    "decreaseInteractiveCount()\0"
};

void SoQtFullViewerP::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        Q_ASSERT(staticMetaObject.cast(_o));
        SoQtFullViewerP *_t = static_cast<SoQtFullViewerP *>(_o);
        switch (_id) {
        case 0: _t->leftWheelPressed(); break;
        case 1: _t->leftWheelChanged((*reinterpret_cast< float(*)>(_a[1]))); break;
        case 2: _t->leftWheelReleased(); break;
        case 3: _t->rightWheelPressed(); break;
        case 4: _t->rightWheelChanged((*reinterpret_cast< float(*)>(_a[1]))); break;
        case 5: _t->rightWheelReleased(); break;
        case 6: _t->bottomWheelPressed(); break;
        case 7: _t->bottomWheelChanged((*reinterpret_cast< float(*)>(_a[1]))); break;
        case 8: _t->bottomWheelReleased(); break;
        case 9: _t->interactbuttonClicked(); break;
        case 10: _t->viewbuttonClicked(); break;
        case 11: _t->homebuttonClicked(); break;
        case 12: _t->sethomebuttonClicked(); break;
        case 13: _t->viewallbuttonClicked(); break;
        case 14: _t->seekbuttonClicked(); break;
        case 15: _t->selectedViewing(); break;
        case 16: _t->selectedDecoration(); break;
        case 17: _t->selectedHeadlight(); break;
        case 18: _t->increaseInteractiveCount(); break;
        case 19: _t->decreaseInteractiveCount(); break;
        default: ;
        }
    }
}

const QMetaObjectExtraData SoQtFullViewerP::staticMetaObjectExtraData = {
    0,  qt_static_metacall 
};

const QMetaObject SoQtFullViewerP::staticMetaObject = {
    { &QObject::staticMetaObject, qt_meta_stringdata_SoQtFullViewerP,
      qt_meta_data_SoQtFullViewerP, &staticMetaObjectExtraData }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &SoQtFullViewerP::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *SoQtFullViewerP::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *SoQtFullViewerP::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_SoQtFullViewerP))
        return static_cast<void*>(const_cast< SoQtFullViewerP*>(this));
    if (!strcmp(_clname, "SoGuiFullViewerP"))
        return static_cast< SoGuiFullViewerP*>(const_cast< SoQtFullViewerP*>(this));
    return QObject::qt_metacast(_clname);
}

int SoQtFullViewerP::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QObject::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 20)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 20;
    }
    return _id;
}
QT_END_MOC_NAMESPACE
