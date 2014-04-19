/****************************************************************************
** Meta object code from reading C++ file 'FacePerformancePlugin.hh'
**
** Created: Fri Apr 18 21:00:18 2014
**      by: The Qt Meta Object Compiler version 63 (Qt 4.8.1)
**
** WARNING! All changes made in this file will be lost!
*****************************************************************************/

#include "../../Plugin-FacePerformance/FacePerformancePlugin.hh"
#if !defined(Q_MOC_OUTPUT_REVISION)
#error "The header file 'FacePerformancePlugin.hh' doesn't include <QObject>."
#elif Q_MOC_OUTPUT_REVISION != 63
#error "This file was generated using the moc from 4.8.1. It"
#error "cannot be used with the include files from this version of Qt."
#error "(The moc has changed too much.)"
#endif

QT_BEGIN_MOC_NAMESPACE
static const uint qt_meta_data_FacePerformancePlugin[] = {

 // content:
       6,       // revision
       0,       // classname
       0,    0, // classinfo
      15,   14, // methods
       0,    0, // properties
       0,    0, // enums/sets
       0,    0, // constructors
       0,       // flags
       7,       // signalCount

 // signals: signature, parameters, type, tag, flags
      23,   22,   22,   22, 0x05,
      54,   36,   22,   22, 0x05,
      98,   83,   22,   22, 0x05,
     128,  119,   22,   22, 0x05,
     155,  141,   22,   22, 0x05,
     204,  184,   22,   22, 0x05,
     246,  232,   22,   22, 0x05,

 // slots: signature, parameters, type, tag, flags
     264,   22,   22,   22, 0x08,
     283,   22,   22,   22, 0x0a,
     301,   22,   22,   22, 0x0a,
     314,   22,   22,   22, 0x0a,
     326,   22,   22,   22, 0x0a,
     345,   22,   22,   22, 0x0a,
     361,   22,   22,   22, 0x0a,
     384,   22,  376,   22, 0x0a,

       0        // eod
};

static const char qt_meta_stringdata_FacePerformancePlugin[] = {
    "FacePerformancePlugin\0\0updateView()\0"
    "_identifier,_type\0updateObject(int,UpdateType)\0"
    "_type,_message\0log(Logtype,QString)\0"
    "_message\0log(QString)\0_name,_widget\0"
    "addToolbox(QString,QWidget*)\0"
    "_filename,_type,_id\0load(QString,DataType,int&)\0"
    "_id,_filename\0save(int,QString)\0"
    "initializePlugin()\0loadBlendshapes()\0"
    "openKinect()\0ReadFrame()\0initRegisterMesh()\0"
    "stopReadFrame()\0trackingMesh()\0QString\0"
    "version()\0"
};

void FacePerformancePlugin::qt_static_metacall(QObject *_o, QMetaObject::Call _c, int _id, void **_a)
{
    if (_c == QMetaObject::InvokeMetaMethod) {
        Q_ASSERT(staticMetaObject.cast(_o));
        FacePerformancePlugin *_t = static_cast<FacePerformancePlugin *>(_o);
        switch (_id) {
        case 0: _t->updateView(); break;
        case 1: _t->updateObject((*reinterpret_cast< int(*)>(_a[1])),(*reinterpret_cast< const UpdateType(*)>(_a[2]))); break;
        case 2: _t->log((*reinterpret_cast< Logtype(*)>(_a[1])),(*reinterpret_cast< QString(*)>(_a[2]))); break;
        case 3: _t->log((*reinterpret_cast< QString(*)>(_a[1]))); break;
        case 4: _t->addToolbox((*reinterpret_cast< QString(*)>(_a[1])),(*reinterpret_cast< QWidget*(*)>(_a[2]))); break;
        case 5: _t->load((*reinterpret_cast< QString(*)>(_a[1])),(*reinterpret_cast< DataType(*)>(_a[2])),(*reinterpret_cast< int(*)>(_a[3]))); break;
        case 6: _t->save((*reinterpret_cast< int(*)>(_a[1])),(*reinterpret_cast< QString(*)>(_a[2]))); break;
        case 7: _t->initializePlugin(); break;
        case 8: _t->loadBlendshapes(); break;
        case 9: _t->openKinect(); break;
        case 10: _t->ReadFrame(); break;
        case 11: _t->initRegisterMesh(); break;
        case 12: _t->stopReadFrame(); break;
        case 13: _t->trackingMesh(); break;
        case 14: { QString _r = _t->version();
            if (_a[0]) *reinterpret_cast< QString*>(_a[0]) = _r; }  break;
        default: ;
        }
    }
}

const QMetaObjectExtraData FacePerformancePlugin::staticMetaObjectExtraData = {
    0,  qt_static_metacall 
};

const QMetaObject FacePerformancePlugin::staticMetaObject = {
    { &QObject::staticMetaObject, qt_meta_stringdata_FacePerformancePlugin,
      qt_meta_data_FacePerformancePlugin, &staticMetaObjectExtraData }
};

#ifdef Q_NO_DATA_RELOCATION
const QMetaObject &FacePerformancePlugin::getStaticMetaObject() { return staticMetaObject; }
#endif //Q_NO_DATA_RELOCATION

const QMetaObject *FacePerformancePlugin::metaObject() const
{
    return QObject::d_ptr->metaObject ? QObject::d_ptr->metaObject : &staticMetaObject;
}

void *FacePerformancePlugin::qt_metacast(const char *_clname)
{
    if (!_clname) return 0;
    if (!strcmp(_clname, qt_meta_stringdata_FacePerformancePlugin))
        return static_cast<void*>(const_cast< FacePerformancePlugin*>(this));
    if (!strcmp(_clname, "BaseInterface"))
        return static_cast< BaseInterface*>(const_cast< FacePerformancePlugin*>(this));
    if (!strcmp(_clname, "ToolboxInterface"))
        return static_cast< ToolboxInterface*>(const_cast< FacePerformancePlugin*>(this));
    if (!strcmp(_clname, "LoggingInterface"))
        return static_cast< LoggingInterface*>(const_cast< FacePerformancePlugin*>(this));
    if (!strcmp(_clname, "LoadSaveInterface"))
        return static_cast< LoadSaveInterface*>(const_cast< FacePerformancePlugin*>(this));
    if (!strcmp(_clname, "OpenFlipper.BaseInterface/1.0"))
        return static_cast< BaseInterface*>(const_cast< FacePerformancePlugin*>(this));
    if (!strcmp(_clname, "OpenFlipper.LoggingInterface/1.0"))
        return static_cast< LoggingInterface*>(const_cast< FacePerformancePlugin*>(this));
    if (!strcmp(_clname, "OpenFlipper.ToolboxInterface/1.1"))
        return static_cast< ToolboxInterface*>(const_cast< FacePerformancePlugin*>(this));
    if (!strcmp(_clname, "OpenFlipper.LoadSaveInterface/1.1"))
        return static_cast< LoadSaveInterface*>(const_cast< FacePerformancePlugin*>(this));
    return QObject::qt_metacast(_clname);
}

int FacePerformancePlugin::qt_metacall(QMetaObject::Call _c, int _id, void **_a)
{
    _id = QObject::qt_metacall(_c, _id, _a);
    if (_id < 0)
        return _id;
    if (_c == QMetaObject::InvokeMetaMethod) {
        if (_id < 15)
            qt_static_metacall(this, _c, _id, _a);
        _id -= 15;
    }
    return _id;
}

// SIGNAL 0
void FacePerformancePlugin::updateView()
{
    QMetaObject::activate(this, &staticMetaObject, 0, 0);
}

// SIGNAL 1
void FacePerformancePlugin::updateObject(int _t1, const UpdateType & _t2)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)), const_cast<void*>(reinterpret_cast<const void*>(&_t2)) };
    QMetaObject::activate(this, &staticMetaObject, 1, _a);
}

// SIGNAL 2
void FacePerformancePlugin::log(Logtype _t1, QString _t2)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)), const_cast<void*>(reinterpret_cast<const void*>(&_t2)) };
    QMetaObject::activate(this, &staticMetaObject, 2, _a);
}

// SIGNAL 3
void FacePerformancePlugin::log(QString _t1)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)) };
    QMetaObject::activate(this, &staticMetaObject, 3, _a);
}

// SIGNAL 4
void FacePerformancePlugin::addToolbox(QString _t1, QWidget * _t2)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)), const_cast<void*>(reinterpret_cast<const void*>(&_t2)) };
    QMetaObject::activate(this, &staticMetaObject, 4, _a);
}

// SIGNAL 5
void FacePerformancePlugin::load(QString _t1, DataType _t2, int & _t3)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)), const_cast<void*>(reinterpret_cast<const void*>(&_t2)), const_cast<void*>(reinterpret_cast<const void*>(&_t3)) };
    QMetaObject::activate(this, &staticMetaObject, 5, _a);
}

// SIGNAL 6
void FacePerformancePlugin::save(int _t1, QString _t2)
{
    void *_a[] = { 0, const_cast<void*>(reinterpret_cast<const void*>(&_t1)), const_cast<void*>(reinterpret_cast<const void*>(&_t2)) };
    QMetaObject::activate(this, &staticMetaObject, 6, _a);
}
QT_END_MOC_NAMESPACE
