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
    using DomainBuilder =  std::function<BaseDomain* (Mesh<Dimension> &, const ElementsIndices &)>;

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

    /**
     * Register an element type to the given VTK cell type and give it a domain builder callback function.
     * @param vtkCellType The type of VTK cell for which we want the builder to be use
     * @param builder     A callback function that will be called with the mesh being created and the list of element
     *                    node indices. The builder is responsible to create the Domain and add it to the mesh. The
     *                    list of element node indices is given as a NxM matrix where N is the number of elements,
     *                    and M is the number of nodes per elements. Domain created domain should be returned by the
     *                    callback, or nullptr if the creation failed.
     *
     *                    DomainBuilder Signature:
     *                    \code{.cpp}
     *                    builder(Mesh<Dimension> & mesh, const ElementsIndices & indices) -> BaseDomain *;
     *                    \endcode
     *
     *                    Where mesh is the current mesh being created and
     *                    ElementsIndices is an array of the following form:
     *
     *                    \verbatim
     *                      // Node 1   Node 2   Node 3    ...    Node M
     *                    [[   e1n1,    e1n2,    e1n3,    ...,    e1nM   ], // Element 1
     *                     [   e2n1,    e2n2,    e2n3,    ...,    e2nM   ], // Element 2
     *                                                    ...
     *                     [   eNn1,    eNn2,    eNn3,    ...,    eNnM]  ] // Element N
     *                     \endverbatim
     *
     *                    Here, for example, e2n4 means the indice of the second node of the fourth element of the mesh.
     *
     * Example:
     * \code{.cpp}
     * auto reader = io::VTKReader<_3D>::Read("meshes/3D_hexahedron_linear.vtk");
     * reader.register_element_type(VTK_HEXAHEDRON, [](Mesh<_3D> & mesh, const Mesh<_3D>::ElementsIndices &indices)) {
     *     // ElementsIndices is a matrix Nx8 : N hexahedral elements of 8 node indices.
     *      return m.add_domain<Hexahedron>("my_hexahdral_domain", indices);
     * });
     * auto mesh = reader.mesh();
     * \endcode
     */
    auto register_element_type(const VTKCellType & vtkCellType, DomainBuilder builder) -> VTKReader & {
        p_domain_builders[vtkCellType] = builder;
        return *this;
    }

private:
    VTKReader(std::string filepath, vtkSmartPointer<vtkUnstructuredGridReader> reader, std::array<UNSIGNED_INTEGER_TYPE, Dimension> axes);

    const std::string p_filepath;
    vtkSmartPointer<vtkUnstructuredGridReader> p_reader;
    std::array<UNSIGNED_INTEGER_TYPE, Dimension> p_axes;
    std::unordered_map<VTKCellType, DomainBuilder> p_domain_builders;
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
