#include <filesystem>
#include <bitset>

#include <Caribou/Constants.h>
#include <Caribou/Geometry/Segment.h>
#include <Caribou/Topology/IO/VTKReader.h>

#include <vtkUnstructuredGrid.h>
#include <vtkPolyData.h>
#include <vtkCellData.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkSmartPointer.h>

namespace fs = std::filesystem;

namespace caribou::topology::io {

template<UNSIGNED_INTEGER_TYPE Dimension>
auto extract_axes_from_3D_vectors(vtkPoints * input_points, const vtkIdType & number_of_points) -> std::array<UNSIGNED_INTEGER_TYPE, Dimension>;

template<UNSIGNED_INTEGER_TYPE Dimension>
VTKReader<Dimension>::VTKReader(std::string filepath, vtkSmartPointer<vtkUnstructuredGridReader> reader, std::array<UNSIGNED_INTEGER_TYPE, Dimension> axes)
: p_filepath(std::move(filepath)), p_reader(std::move(reader)), p_axes(axes)
{
    // Segments
    register_element_type<geometry::Segment<Dimension, Linear>>(VTK_LINE);
    register_element_type<geometry::Segment<Dimension, Quadratic>>(VTK_QUADRATIC_EDGE);
}

template<UNSIGNED_INTEGER_TYPE Dimension>
auto VTKReader<Dimension>::Read(const std::string &filepath) -> VTKReader<Dimension> {
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


    return VTKReader<Dimension>(filepath, reader, axes);
}

template<UNSIGNED_INTEGER_TYPE Dimension>
auto VTKReader<Dimension>::mesh () const -> UnstructuredMesh<Dimension> {
    using WorldCoordinates = typename UnstructuredMesh<Dimension>::WorldCoordinates;

    UnstructuredMesh<Dimension> m;
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

    m = UnstructuredMesh<Dimension> (nodes);

    // Import elements
    vtkSmartPointer <vtkCellTypes> types = vtkSmartPointer <vtkCellTypes>::New();
    p_reader->GetOutput()->GetCellTypes(types);
    vtkIdType number_of_element_types = types->GetNumberOfTypes();
    for (unsigned int i = 0; i < number_of_element_types; ++i) {
        const auto type = types->GetCellType(i);

        if (p_domain_builders.find(static_cast<VTKCellType>(type)) == p_domain_builders.end()) {
            // This element type isn't supported (no domain builder found)
            continue;
        }

        const auto & domain_builder = p_domain_builders.at(static_cast<VTKCellType>(type));
        BaseDomain * domain = domain_builder(m);

        auto cells = vtkSmartPointer <vtkIdTypeArray>::New();
        p_reader->GetOutput()->GetIdsOfCellsOfType(type, cells);
        for (vtkIdType j = 0; j < cells->GetDataSize(); ++j) {
            vtkIdType cell_id = cells->GetValue(j);
            vtkCell* cell = p_reader->GetOutput()->GetCell(cell_id);
            std::vector<UNSIGNED_INTEGER_TYPE> node_indices (cell->GetNumberOfPoints());
            for (vtkIdType k = 0; k < cell->GetNumberOfPoints(); ++k) {
                node_indices[k] = static_cast<UNSIGNED_INTEGER_TYPE>(cell->GetPointId(k));
            }
            domain->add_element(node_indices.data(), node_indices.size());
        }
    }


    return m;
}


template<UNSIGNED_INTEGER_TYPE Dimension>
void VTKReader<Dimension>::print (std::ostream & out) {
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
 *    [[24, 0, 12], [64, 0, 22], [51, 0, 9]]
 * will produce the following axis for a 2D mesh
 *    [0, 2] // x is the 0-axis, y is the 2-axis
 * and will raise an error for a 1D mesh
 *
 * The input
 *    [[5, 0, 12], [5, 0, 22], [5, 0, 9]]
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

        for (std::size_t axis = 0, c = 0; axis < 3; ++axis) {
            if (not is_all_the_same[axis]) {
                axes[c++] = axis;
            }
        }
    }

    return axes;
}

template class VTKReader<1>;
template class VTKReader<2>;
template class VTKReader<3>;
}