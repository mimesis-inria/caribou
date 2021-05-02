#include <filesystem>
#include <bitset>

#include <Caribou/constants.h>
#include <Caribou/Geometry/Quad.h>
#include <Caribou/Geometry/Triangle.h>
#include <Caribou/Geometry/Segment.h>
#include <Caribou/Geometry/Tetrahedron.h>
#include <Caribou/Geometry/Hexahedron.h>

#include <Caribou/Topology/IO/VTKReader.h>

#include <vtkUnstructuredGrid.h>
#include <vtkPolyData.h>
#include <vtkCellData.h>
#include <vtkCellType.h>
#include <vtkCellTypes.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkSmartPointer.h>
#include <vtkPointData.h>

namespace fs = std::filesystem;

namespace caribou::topology::io {

template<UNSIGNED_INTEGER_TYPE Dimension>
auto extract_axes_from_3D_vectors(vtkPoints * input_points, const vtkIdType & number_of_points) -> std::array<UNSIGNED_INTEGER_TYPE, Dimension>;

template<UNSIGNED_INTEGER_TYPE Dimension, typename NodeIndex>
VTKReader<Dimension, NodeIndex>::VTKReader(std::string filepath, vtkSmartPointer<vtkUnstructuredGridReader> reader, std::array<UNSIGNED_INTEGER_TYPE, Dimension> axes)
: p_filepath(std::move(filepath)), p_reader(std::move(reader)), p_axes(axes)
{
    // Segments
    register_element_type<geometry::Segment<Dimension, Linear>>(VTK_LINE);
    register_element_type<geometry::Segment<Dimension, Quadratic>>(VTK_QUADRATIC_EDGE);

    if constexpr (Dimension > 1) {
        // Quads
        register_element_type<geometry::Quad<Dimension, Linear>>(VTK_QUAD);
        register_element_type<geometry::Quad<Dimension, Quadratic>>(VTK_QUADRATIC_QUAD);

        // Triangles
        register_element_type<geometry::Triangle<Dimension, Linear>>(VTK_TRIANGLE);
        register_element_type<geometry::Triangle<Dimension, Quadratic>>(VTK_QUADRATIC_TRIANGLE);
    }

    if constexpr (Dimension > 2) {
        // Tetrahedrons
        register_element_type<geometry::Tetrahedron<Linear>>(VTK_TETRA);
        register_element_type<geometry::Tetrahedron<Quadratic>>(VTK_QUADRATIC_TETRA);

        // Hexahedrons
        register_element_type<geometry::Hexahedron<Linear>>(VTK_HEXAHEDRON);
        register_element_type(VTK_QUADRATIC_HEXAHEDRON, [](Mesh<Dimension> & m, const ElementsIndices & indices) {

            // The order of quadratic hexahedron node indces in VTK aren't the same as the order used in Caribou
            ElementsIndices reordered_indices = indices;
//            reordered_indices.col(9)  = indices.col(9);
            reordered_indices.col(10) = indices.col(18);
//            reordered_indices.col(11) = indices.col(11);
            reordered_indices.col(12) = indices.col(10);
            reordered_indices.col(13) = indices.col(14);
            reordered_indices.col(14) = indices.col(15);
            reordered_indices.col(15) = indices.col(12);
            reordered_indices.col(16) = indices.col(13);
            reordered_indices.col(17) = indices.col(16);
            reordered_indices.col(18) = indices.col(19);
            reordered_indices.col(19) = indices.col(17);
            return m.template add_domain<geometry::Hexahedron<Quadratic>, NodeIndex>("domain_"+std::to_string(m.number_of_domains()+1), reordered_indices);
        });
    }
}

template<UNSIGNED_INTEGER_TYPE Dimension, typename NodeIndex>
auto VTKReader<Dimension, NodeIndex>::Read(const std::string &filepath) -> VTKReader<Dimension, NodeIndex> {
    // Make sure filepath is a valid path
    if (not fs::exists(filepath)) {
        throw std::invalid_argument("File '" + filepath + "' does not exists or cannot be read.");
    }

    // Get all data from the file
    auto reader = vtkSmartPointer<vtkUnstructuredGridReader>::New();
    reader->SetFileName(filepath.c_str());
    reader->Update();

    vtkUnstructuredGrid * output = reader->GetOutput();

    const auto axes = extract_axes_from_3D_vectors<Dimension>(output->GetPoints(), output->GetNumberOfPoints());


    return VTKReader<Dimension, NodeIndex>(filepath, reader, axes);
}

template<UNSIGNED_INTEGER_TYPE Dimension, typename NodeIndex>
auto VTKReader<Dimension, NodeIndex>::mesh () const -> Mesh<Dimension> {
    using WorldCoordinates = typename Mesh<Dimension>::WorldCoordinates;

    Mesh<Dimension> m;
    const auto number_of_nodes = p_reader->GetOutput()->GetNumberOfPoints();
    if (number_of_nodes == 0) {
        return m;
    }

    // Import nodes
    std::vector<WorldCoordinates> nodes;
    nodes.resize(p_reader->GetOutput()->GetNumberOfPoints());

    for (vtkIdType i = 0; i < number_of_nodes; ++i) {
        double * v = p_reader->GetOutput()->GetPoints()->GetPoint(i);
        for (std::size_t axis = 0; axis < Dimension; ++axis) {
            nodes[i][axis] = v[p_axes[axis]];
        }
    }

    m = Mesh<Dimension> (nodes);

    // Import elements
    vtkSmartPointer <vtkCellTypes> types = vtkSmartPointer <vtkCellTypes>::New();
    p_reader->GetOutput()->GetCellTypes(types);
    vtkIdType number_of_element_types = types->GetNumberOfTypes();
    for (unsigned int i = 0; i < number_of_element_types; ++i) {
        const auto type = types->GetCellType(i);
        auto cells = vtkSmartPointer <vtkIdTypeArray>::New();
        p_reader->GetOutput()->GetIdsOfCellsOfType(type, cells);
        const auto number_of_elements = cells->GetDataSize();

        if (number_of_elements == 0) {
            continue;
        }

        if (p_domain_builders.find(static_cast<VTKCellType>(type)) == p_domain_builders.end()) {
            // This element type isn't supported (no domain builder found)
            continue;
        }

        const auto number_of_nodes_per_element = p_reader->GetOutput()->GetCell(cells->GetValue(0))->GetNumberOfPoints();

        Eigen::Matrix<NodeIndex, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> indices;
        indices.resize(number_of_elements, number_of_nodes_per_element);

        for (vtkIdType j = 0; j < number_of_elements; ++j) {
            vtkIdType cell_id = cells->GetValue(j);
            vtkCell* cell = p_reader->GetOutput()->GetCell(cell_id);
            std::vector<NodeIndex> node_indices (cell->GetNumberOfPoints());
            for (vtkIdType k = 0; k < cell->GetNumberOfPoints(); ++k) {
                indices(j,k) = static_cast<NodeIndex>(cell->GetPointId(k));
            }
        }

        const auto & domain_builder = p_domain_builders.at(static_cast<VTKCellType>(type));
        domain_builder(m, indices);
    }

    return m;
}


template<UNSIGNED_INTEGER_TYPE Dimension, typename NodeIndex>
void VTKReader<Dimension, NodeIndex>::print (std::ostream & out) const {
    vtkUnstructuredGrid * output = p_reader->GetOutput();

    out << "input has " << output->GetNumberOfPoints() << " points.\n";
    out << "input has " << output->GetNumberOfCells() << " cells.\n";

    vtkSmartPointer <vtkCellTypes> types = vtkSmartPointer <vtkCellTypes>::New();
    output->GetCellTypes(types);
    out << types->GetNumberOfTypes() << " types\n";
    for (unsigned int i = 0; i < types->GetNumberOfTypes(); ++i) {
        const auto type = types->GetCellType(i);
        auto cells = vtkSmartPointer <vtkIdTypeArray>::New();
        output->GetIdsOfCellsOfType(type, cells);
        out << vtkCellTypes::GetClassNameFromTypeId(type) << " : " << cells->GetDataSize() << "\n";
    }

    vtkIdType numberOfCellArrays = output->GetCellData()->GetNumberOfArrays();
    if (numberOfCellArrays > 0) {
        out << "Cell data arrays:\n";
        for (vtkIdType i = 0; i < numberOfCellArrays; i++) {
            const auto dataTypeID = output->GetCellData()->GetArray(i)->GetDataType();
            const auto dataTypeIDStr = output->GetCellData()->GetArray(i)->GetDataTypeAsString();
            out << "  Array " << i << ": "
                << output->GetCellData()->GetArrayName(i)
                << "  (type: " << dataTypeIDStr << " - " << dataTypeID << ")\n";
        }
    }

    vtkIdType numberOfFieldArrays = output->GetFieldData()->GetNumberOfArrays();
    if (numberOfFieldArrays) {
        out << "Field data arrays:\n";
        for (vtkIdType i = 0; i < numberOfFieldArrays; i++) {
            const auto dataTypeID = output->GetFieldData()->GetArray(i)->GetDataType();
            const auto dataTypeIDStr = output->GetFieldData()->GetArray(i)->GetDataTypeAsString();
            out << "  Array " << i << ": "
                << output->GetFieldData()->GetArrayName(i)
                << "  (type: " << dataTypeIDStr << " - " << dataTypeID << ")\n";
        }
    }
}

/**
 * Extract the axes from an array of 3D coordinates.
 *
 * This method will first check which axis (x, y or z) has all of its coordinates equals.
 *
 * For example, the input
 *    [[24, 0, 12],
 *     [64, 0, 22],
 *     [51, 0, 9]]
 * will produce the following axis for a 2D mesh
 *    [0, 2] // x is the 0-axis, y is the 2-axis
 * and will raise an error for a 1D mesh
 *
 * The input
 *    [[5, 0, 12],
 *     [5, 0, 22],
 *     [5, 0, 9]]
 * will produce the following axis for a 1D mesh
 *    [2] // x is the 2-axis
 * and will raise an error for a 2D mesh
 *
 * @tparam Dimension The dimension of the output coordinates]
 * @param input_points Input array of 3D vectors
 * @param number_of_points Number of points in the input array
 */
template<UNSIGNED_INTEGER_TYPE Dimension>
auto extract_axes_from_3D_vectors(vtkPoints * input_points, const vtkIdType & number_of_points) -> std::array<UNSIGNED_INTEGER_TYPE, Dimension>
{
    // Find the good axes (only for dimensions 1 and 2)
    std::array<UNSIGNED_INTEGER_TYPE, Dimension> axes {};
    if (Dimension == 3) {
      axes[0] = 0;
      axes[1] = 1;
      axes[2] = 2;
    } else {
        double * v = input_points->GetPoint(0);
        std::bitset<3> is_all_the_same(0b111); // Set all to '1' at first
        std::array<double, 3> last_value {v[0], v[1], v[2]};
        for (vtkIdType i = 1; i < number_of_points; ++i) {
            v = input_points->GetPoint(i);
            for (std::size_t axis = 0; axis < 3; ++axis) {
                if (last_value[axis] != v[axis])
                    is_all_the_same[axis] = false;
            }
        }

        if (is_all_the_same.flip().count() != Dimension) {
            throw std::runtime_error("Unable to convert a 3D field to a " + std::to_string(Dimension) +
                                     "D field. Their is " + std::to_string(is_all_the_same.count()) +
                                     " axes from the input mesh that have the same value.");
        }

        // Flip back
        is_all_the_same.flip();

        for (std::size_t axis = 0, c = 0; axis < 3; ++axis) {
            if (not is_all_the_same[axis]) {
                axes[c++] = axis;
            }
        }
    }

    return axes;
}

template class VTKReader<1, unsigned long long int>;
template class VTKReader<1, unsigned long int>;
template class VTKReader<1, unsigned int>;

template class VTKReader<2, unsigned long long int>;
template class VTKReader<2, unsigned long int>;
template class VTKReader<2, unsigned int>;

template class VTKReader<3, unsigned long long int>;
template class VTKReader<3, unsigned long int>;
template class VTKReader<3, unsigned int>;
}