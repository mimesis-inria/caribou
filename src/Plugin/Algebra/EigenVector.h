#pragma once

#include <Caribou/config.h>
#include <Caribou/macros.h>
#include <Caribou/traits.h>
#include <sofa/defaulttype/BaseMatrix.h>

#include <Eigen/Dense>

namespace SofaCaribou::Algebra {

/**
 * This class can be use to represent any Eigen vector within SOFA. It
 * implements (and overloads) the BaseVector class of Sofa. This class hence
 * allows to either to create a new Eigen vector or copy/map an existing one
 * and use it in Sofa's components or visitors.
 *
 * @tparam Derived
 * @tparam Enable
 */
template <typename Derived, typename Enable = void>
class EigenVector : public sofa::defaulttype::BaseVector
{
    static_assert(
            std::is_base_of_v<Eigen::EigenBase<std::decay_t<Derived> >, std::decay_t<Derived> >,
            "The class template argument 'Derived' must be derived from Eigen::EigenBase<Derived>."
    );

    static_assert(
        std::remove_reference_t<Derived>::IsVectorAtCompileTime,
        "The class template argument 'Derived' must be an Eigen vector (either 1 or dynamic column"
        " at compile time, or 1 or dynamic row at compile time."
    );

public:
    using EigenType = std::remove_cv_t<std::remove_reference_t<Derived>>;
    using Base = sofa::defaulttype::BaseVector;
    using Index = Base::Index;
    using Real = typename EigenType::Scalar;

    /**
    * Construct the class using another Eigen vector. Depending on the template parameter used for the class,
    * this constructor will either create a new Eigen vector and copy its content from the parameter eigen_matrix,
    * or it will simply store a reference to this external matrix.
    * @param eigen_matrix The external matrix
    */
    explicit EigenVector(std::remove_reference_t<Derived> & eigen_vector) : p_eigen_vector(eigen_vector) {}

    /**
     * Construct a new Eigen vector of type Derived have n elements.
     * @param n Number of elements in the vector.
     */
    EigenVector(Eigen::Index n) : p_eigen_vector(n) {}

    /** Number of elements in the vector */
    [[nodiscard]]
    Index size() const final {return static_cast<int>(p_eigen_vector.size());}

    /** Read the value of element i */
    [[nodiscard]]
    SReal element(Index index) const final {return static_cast<SReal>(p_eigen_vector[index]);}

    /** Resize the vector. This method resets to zero all entries. */
    void resize(Index n) final {p_eigen_vector.resize(n);p_eigen_vector.setZero();}

    /** Reset all values to 0 */
    void clear() final {p_eigen_vector.setZero();}

    /** Write the value of element index */
    void set(Index index, SReal value) final {p_eigen_vector[index] = static_cast<Real>(value);}

    /** Add v to the existing value of element i */
    void add(Index index, SReal value) final {p_eigen_vector[index] += static_cast<Real>(value);}

    /** Type of elements stored in this matrix */
    Base::ElementType getElementType() const final {
        return (std::is_integral_v<Real> ? Base::ELEMENT_INT : Base::ELEMENT_FLOAT);
    }

    /** Size of one element stored in this matrix.*/
    std::size_t getElementSize() const final {
        return sizeof(Real);
    }

    /** Get a const reference to the underlying Eigen vector  */
    const EigenType & vector() const {return p_eigen_vector;}

    /** Get a reference to the underlying Eigen vector  */
    EigenType & vector() {return p_eigen_vector;}

private:
    ///< The actual Eigen Vector.
    Derived p_eigen_vector;
};

} // namespace SofaCaribou::Algebra
