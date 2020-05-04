

#include <pybind11/pybind11.h>
#include <lm/pylm.h>
#include "arepoloader.h"

LM_NAMESPACE_BEGIN(LM_NAMESPACE)
/*


PYBIND11_MODULE(arepopy, m) {
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
        

        virtual Bound bound() const override  {
            PYBIND11_OVERLOAD(Bound, Volume_Arepo, bound);
        }
        virtual Float max_scalar() const override {
            PYBIND11_OVERLOAD(Float, Volume_Arepo, max_scalar);
        }

        virtual bool has_scalar() const override {
            PYBIND11_OVERLOAD(bool, Volume_Arepo, has_scalar);
        }

        virtual Float eval_scalar(Vec3 p) const override {
            PYBIND11_OVERLOAD(Float, Volume_Arepo, eval_scalar, p);
        }

        virtual Vec3 eval_color(Vec3 p) const override {
            PYBIND11_OVERLOAD(Vec3, Volume_Arepo, eval_color, p);
        }

        virtual bool has_color() const override {
            PYBIND11_OVERLOAD(bool, Volume_Arepo, has_color);
        }

        //TODO other functions ?
    };
    std::cout<<"executed"<<std::endl;
    // You must not add .def() for construct() function
    // which is already registered in parent class.
    pybind11::class_<Volume_Arepo, VolumeArepo_Py, Component::Ptr<Volume_Arepo>>(m, "Volume_Arepo")
        .def(pybind11::init<>())
        .def("has_color",&Volume_Arepo::has_color)
        .def("eval_color", &Volume_Arepo::eval_color)
        .def("eval_scalar", &Volume_Arepo::eval_scalar)
        .def("has_scalar",&Volume_Arepo::has_scalar)
        .def("max_scalar", &Volume_Arepo::max_scalar)
        .def("bound", &Volume_Arepo::bound)
        //.def("feedback", &Tuner::feedback)
        //.def("getConf", &Tuner::getConf)
        .PYLM_DEF_COMP_BIND(Volume_Arepo);
    std::cout << "executed" << std::endl;
}

*/

LM_NAMESPACE_END(LM_NAMESPACE)
