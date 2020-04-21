#pragma once

#include <string>
#include <iostream>
#include <utility>
#include <functional>
#include <unordered_map>

#include <Caribou/config.h>
#include <Caribou/Topology/Mesh.h>
#include <Caribou/Topology/Domain.h>

#include <vtkSmartPointer.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkCellType.h>

namespace caribou::topology::io {

template<UNSIGNED_INTEGER_TYPE Dimension>
class VTKReader {
    static_assert(Dimension == 1 or Dimension == 2 or Dimension == 3, "The VTKReader can only read 1D, 2D or 3D fields.");
public:

    using ElementsIndices = Eigen::Matrix<UNSIGNED_INTEGER_TYPE, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

    /** Build a new VTKReader instance by reading a vtk file. */
    static auto Read(const std::string & filepath) -> VTKReader;

    /** Print information about the current vtk file. */
    void print (std::ostream &out) const;

    /** Build the actual unstructured mesh from the vtk file. */
    [[nodiscard]]
    auto mesh() const -> Mesh<Dimension>;

    /** Register an element type to the given VTK cell type */
    template<typename Element>
    auto register_element_type(const VTKCellType & vtkCellType) -> VTKReader & {
        p_domain_builders[vtkCellType] = [](Mesh<Dimension> & m, const ElementsIndices & indices) {
            return m.template add_domain<Element>("domain_"+std::to_string(m.number_of_domains()+1), indices);
        };
        return *this;
    }

private:
    VTKReader(std::string filepath, vtkSmartPointer<vtkUnstructuredGridReader> reader, std::array<UNSIGNED_INTEGER_TYPE, Dimension> axes);

    const std::string p_filepath;
    vtkSmartPointer<vtkUnstructuredGridReader> p_reader;
    std::array<UNSIGNED_INTEGER_TYPE, Dimension> p_axes;
    std::unordered_map<VTKCellType, std::function<BaseDomain* (Mesh<Dimension> &, const ElementsIndices &)>> p_domain_builders;
};

extern template class VTKReader<1>;
extern template class VTKReader<2>;
extern template class VTKReader<3>;

template<UNSIGNED_INTEGER_TYPE Dimension>
auto operator<<(std::ostream& os, const caribou::topology::io::VTKReader<Dimension> & t) -> std::ostream&
{
    t.print(os);
    return os;
}

} /// namespace caribou::topology::io
