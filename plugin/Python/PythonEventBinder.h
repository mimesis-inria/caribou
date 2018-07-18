#ifndef SOFAMESHLESS_PYTHONEVENTBINDER_H
#define SOFAMESHLESS_PYTHONEVENTBINDER_H

#include <sofa/core/objectmodel/BaseObject.h>
#include <sofa/core/objectmodel/Event.h>
#include <SofaPython/PythonCommon.h>

namespace sofa {
namespace caribou {
namespace python {

class PythonEventBinder : public virtual sofa::core::objectmodel::BaseObject
{
public:
    SOFA_CLASS(PythonEventBinder, sofa::core::objectmodel::BaseObject);

    PythonEventBinder() {
        this->f_listening.setValue(true);
    };

    void handleEvent(sofa::core::objectmodel::Event* event) override;
    void bind(std::string event_name, PyObject* o);
    void bind(size_t event_index, PyObject* o);
    void raise(size_t event_index);

private:
    sofa::helper::vector<std::vector<PyObject*>> m_callbacks;
};

}
}
}

#endif //SOFAMESHLESS_PYTHONEVENTBINDER_H
