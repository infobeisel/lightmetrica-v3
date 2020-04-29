#include "arepoloader.h"


#include <pybind11/pybind11.h>
#include <lm/pylm.h>

LM_NAMESPACE_BEGIN(LM_NAMESPACE)



PYBIND11_MODULE(py_arepoloader, m) {
    m.doc() = "Python binding illustris arepo";
    std::cout<<"executed"<<std::endl;
    // Bind free function
    //m.def("test_function", &test_function);
    
    // Define a trampoline class for the interface
    // ref. https://pybind11.readthedocs.io/en/stable/advanced/classes.html
    class VolumeArepo_Py final : public Volume_Arepo {
        virtual void construct(const lm::Json& prop) override {
            LM_INFO("try to construct Volume Arepo");

            PYBIND11_OVERLOAD(void, Volume_Arepo, construct, prop);
        }
        //TODO other functions ?
    };
    std::cout<<"executed"<<std::endl;
    // You must not add .def() for construct() function
    // which is already registered in parent class.
    pybind11::class_<Volume_Arepo, VolumeArepo_Py, Component::Ptr<Volume_Arepo>>(m, "Volume_Arepo")
        .def(pybind11::init<>())
        //.def("feedback", &Tuner::feedback)
        //.def("getConf", &Tuner::getConf)
        .PYLM_DEF_COMP_BIND(Volume_Arepo);
    std::cout << "executed" << std::endl;
}



LM_NAMESPACE_END(LM_NAMESPACE)
