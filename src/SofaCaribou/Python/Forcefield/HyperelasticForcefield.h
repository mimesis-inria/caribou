#pragma once

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>

#include <SofaPython3/PythonFactory.h>
#include <SofaPython3/Sofa/Core/Binding_Base.h>
#include <SofaPython3/Sofa/Core/Data/Binding_DataContainer.h>
#include <SofaPython3/PythonEnvironment.h>

#include <SofaCaribou/Forcefield/HyperelasticForcefield.h>
#include <SofaCaribou/Forcefield/HyperelasticForcefield.inl>

#include <sofa/core/MechanicalParams.h>
#include <sofa/core/behavior/MultiMatrixAccessor.h>

using namespace pybind11::literals;

namespace SofaCaribou::forcefield::python {

template <typename Element>
class HyperelasticForcefieldTrampoline : public HyperelasticForcefield<Element> {
public:
    SOFA_CLASS(SOFA_TEMPLATE(HyperelasticForcefieldTrampoline, Element), SOFA_TEMPLATE(HyperelasticForcefield, Element));
    using typename HyperelasticForcefield<Element>::DataVecDeriv;
    using typename HyperelasticForcefield<Element>::DataVecCoord;

    void init() override {
        sofapython3::PythonEnvironment::gil acquire;
        PYBIND11_OVERLOAD(void, HyperelasticForcefield<Element>, init);
    }

    void addForce(const sofa::core::MechanicalParams* mparams,  DataVecDeriv& f, const DataVecCoord& x, const DataVecDeriv& v) override {
        sofapython3::PythonEnvironment::gil acquire;

        // Try to look up the overridden method on the Python side.
#if PYBIND11_VERSION_MAJOR <= 2 && PYBIND11_VERSION_MINOR < 6
        pybind11::function override = pybind11::get_overload(this, "addForce");
#else
        pybind11::function override = pybind11::get_override(this, "addForce");
#endif
        if (override) {
            // python override method is found, call it
            pybind11::dict mp = pybind11::dict("time"_a=this->getContext()->getTime(),
                                               "mFactor"_a=mparams->mFactor(),
                                               "bFactor"_a=mparams->bFactor(),
                                               "kFactor"_a=mparams->kFactor(),
                                               "isImplicit"_a=mparams->implicit(),
                                               "energy"_a=mparams->energy());
            override(mp, sofapython3::PythonFactory::toPython(&f), sofapython3::PythonFactory::toPython(&x), sofapython3::PythonFactory::toPython(&v));
        } else {
            // No python override, call the generic method
            HyperelasticForcefield<Element>::addForce(mparams, f, x, v);
        }
    }

    void addDForce(const sofa::core::MechanicalParams * mparams, DataVecDeriv & df, const DataVecDeriv & dx) override {
        sofapython3::PythonEnvironment::gil acquire;

        // Try to look up the overridden method on the Python side.
#if PYBIND11_VERSION_MAJOR <= 2 && PYBIND11_VERSION_MINOR < 6
        pybind11::function override = pybind11::get_overload(this, "addDForce");
#else
        pybind11::function override = pybind11::get_override(this, "addDForce");
#endif
        if (override) {
            // python override method is found, call it
            pybind11::dict mp = pybind11::dict("time"_a=this->getContext()->getTime(),
                                               "mFactor"_a=mparams->mFactor(),
                                               "bFactor"_a=mparams->bFactor(),
                                               "kFactor"_a=mparams->kFactor(),
                                               "isImplicit"_a=mparams->implicit(),
                                               "energy"_a=mparams->energy());
            override(mp, sofapython3::PythonFactory::toPython(&df), sofapython3::PythonFactory::toPython(&dx));
        } else {
            // No python override, call the generic method
            HyperelasticForcefield<Element>::addDForce(mparams, df, dx);
        }
    }
};

template<typename Element>
void bind_hyperelastic_forcefield(pybind11::module &m, const std::string & template_name) {
    pybind11::module::import("Sofa");

    std::string name = "HyperelasticForcefield_" + template_name;

    using Real = typename HyperelasticForcefield<Element>::Real;
    using VecCoord = typename HyperelasticForcefield<Element>::VecCoord;
    using VecDeriv = typename HyperelasticForcefield<Element>::VecDeriv;


    pybind11::class_<HyperelasticForcefield<Element>, sofa::core::objectmodel::BaseObject, HyperelasticForcefieldTrampoline<Element>, sofapython3::py_shared_ptr<HyperelasticForcefield<Element>>> c(m, name.c_str());
    c.def("init", &HyperelasticForcefield<Element>::init);
    c.def("K", &HyperelasticForcefield<Element>::K);
    c.def("cond", &HyperelasticForcefield<Element>::cond);
    c.def("eigenvalues", &HyperelasticForcefield<Element>::eigenvalues);
    c.def("assemble_stiffness", [](HyperelasticForcefield<Element> & self, const Eigen::Matrix<double, Eigen::Dynamic, HyperelasticForcefield<Element>::Dimension, Eigen::RowMajor> & x) {
        self.assemble_stiffness(x);
    }, pybind11::arg("x").noconvert(true));
    c.def("assemble_stiffness", [](HyperelasticForcefield<Element> & self, const Eigen::Matrix<float, Eigen::Dynamic, HyperelasticForcefield<Element>::Dimension, Eigen::RowMajor> & x) {
        self.assemble_stiffness(x);
    }, pybind11::arg("x").noconvert(true));
    c.def("addForce",
          [](HyperelasticForcefield<Element> &self,
             const pybind11::dict &mp,
             sofapython3::DataContainer &f,
             const sofapython3::DataContainer &x,
             const sofapython3::DataContainer &v
          ) {
              auto *ff = reinterpret_cast<sofa::core::objectmodel::Data<VecDeriv> *>(&f);
              const auto *xx = reinterpret_cast<const sofa::core::objectmodel::Data<VecCoord> *>(&x);
              const auto *vv = reinterpret_cast<const sofa::core::objectmodel::Data<VecDeriv> *>(&v);
              sofa::core::MechanicalParams mparams;
              if (mp.template contains("mFactor")) {
                  mparams.setMFactor(pybind11::cast<Real>(mp["mFactor"]));
              }
              if (mp.template contains("bFactor")) {
                  mparams.setBFactor(pybind11::cast<Real>(mp["bFactor"]));
              }
              if (mp.template contains("kFactor")) {
                  mparams.setKFactor(pybind11::cast<Real>(mp["kFactor"]));
              }
              if (mp.template contains("isImplicit")) {
                  mparams.setImplicit(pybind11::cast<bool>(mp["isImplicit"]));
              }
              if (mp.template contains("energy")) {
                  mparams.setEnergy(pybind11::cast<bool>(mp["energy"]));
              }
              if (mp.template contains("dt")) {
                  mparams.setDt(pybind11::cast<Real>(mp["dt"]));
              }

              self.addForce(&mparams, *ff, *xx, *vv);
          },
          pybind11::arg("mechanical_params"),
          pybind11::arg("f").noconvert(),
          pybind11::arg("x").noconvert(),
          pybind11::arg("v").noconvert()
    );

    c.def("addDForce",
          [](HyperelasticForcefield<Element> & self,
                  const pybind11::dict & mp,
                  sofapython3::DataContainer &df,
                  const sofapython3::DataContainer &dx
                  )
                  {
        auto * dff = reinterpret_cast<sofa::core::objectmodel::Data<VecDeriv>*>(&df);
        const auto * dxx = reinterpret_cast<const sofa::core::objectmodel::Data<VecCoord>*>(&dx);
        sofa::core::MechanicalParams mparams;
        if (mp.template contains("mFactor")) {
            mparams.setMFactor(pybind11::cast<Real>(mp["mFactor"]));
        }
        if (mp.template contains("bFactor")) {
            mparams.setBFactor(pybind11::cast<Real>(mp["bFactor"]));
        }
        if (mp.template contains("kFactor")) {
            mparams.setKFactor(pybind11::cast<Real>(mp["kFactor"]));
        }
        if (mp.template contains("isImplicit")) {
            mparams.setImplicit(pybind11::cast<bool>(mp["isImplicit"]));
        }
        if (mp.template contains("energy")) {
            mparams.setEnergy(pybind11::cast<bool>(mp["energy"]));
        }
        if (mp.template contains("dt")) {
            mparams.setDt(pybind11::cast<Real>(mp["dt"]));
        }

        self.addDForce(&mparams, *dff, *dxx);
        },
        pybind11::arg("mechanical_params"),
        pybind11::arg("df").noconvert(),
        pybind11::arg("dx").noconvert()
        );

    c.def(pybind11::init([](pybind11::args &args, pybind11::kwargs &kwargs) {
        auto ff = sofa::core::sptr<HyperelasticForcefieldTrampoline<Element>> (new HyperelasticForcefieldTrampoline<Element>());

        ff->f_listening.setValue(true);

        if (args.size() == 1) ff->setName(pybind11::cast<std::string>(args[0]));

        pybind11::object cc = pybind11::cast(ff);
        for (auto kv : kwargs) {
            auto key = pybind11::cast<std::string>(kv.first);
            auto value = pybind11::reinterpret_borrow<pybind11::object>(kv.second);
            if (key == "name") {
                if (!args.empty()) {
                    throw pybind11::type_error("The name is set twice as a "
                                               "named argument='" + pybind11::cast<std::string>(value) + "' and as a"
                                                                                             "positional argument='" +
                                                                                             pybind11::cast<std::string>(args[0]) + "'.");
                }
            }
            sofapython3::BindingBase::SetAttr(cc, key, value);
        }
        return ff;
    }));

    sofapython3::PythonFactory::registerType<HyperelasticForcefield<Element>>([template_name](sofa::core::objectmodel::Base* o) {
        return pybind11::cast(dynamic_cast<HyperelasticForcefield<Element>*>(o));
    });
}

template<typename Element>
void bind_hyperelastic_forcefield(pybind11::module &m) {
    bind_hyperelastic_forcefield<Element>(m, HyperelasticForcefield<Element>::templateName());
}

void addHyperElasticForcefield(pybind11::module &m);
}