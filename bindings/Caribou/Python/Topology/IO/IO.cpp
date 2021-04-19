#include <pybind11/pybind11.h>

#include <Caribou/Topology/config.h>

#ifdef CARIBOU_WITH_VTK
#include <Caribou/Topology/IO/VTKReader.h>
#endif
namespace caribou::topology::io::bindings {

#ifdef CARIBOU_WITH_VTK
template<UNSIGNED_INTEGER_TYPE Dimension>
void add_reader(pybind11::module & m) {
    std::string name = "VTKReader" + std::to_string(Dimension) + "D";
    pybind11::class_<VTKReader<Dimension>> c(m, name.c_str());

    c.def_static("Read", &VTKReader<Dimension>::Read);
    c.def("__str__", [](const VTKReader<Dimension> & self) {
        std::stringstream ss;
        ss << self;
        return ss.str();
    });
}
#endif

void create_IO(pybind11::module & m) {
    pybind11::module io = m.def_submodule("IO");

#ifdef CARIBOU_WITH_VTK
    add_reader<1>(m);
    add_reader<2>(m);
    add_reader<3>(m);

    io.def("VTKReader", [](const std::string & filepath, unsigned int dimension) {
        if (dimension == 1) {
            return pybind11::cast(VTKReader<1>::Read(filepath));
        } else if (dimension == 2) {
            return pybind11::cast(VTKReader<2>::Read(filepath));
        } else if (dimension == 3) {
            return pybind11::cast(VTKReader<3>::Read(filepath));
        } else {
            throw std::runtime_error("Trying to create a VTKReader with a dimension that is not 1, 2 or 3.");
        }
    }, pybind11::arg("filepath"), pybind11::arg("dimension") = 3);
#endif
}

} // namespace caribou::topology::io::bindings