#ifndef CARIBOU_ALGEBRA_VECTOR_H
#define CARIBOU_ALGEBRA_VECTOR_H

#include <cstddef>
#include <array>
#include <initializer_list>
#include <algorithm>
#include <numeric>
#include <cmath>

namespace caribou
{
namespace algebra
{

/**
 * A simple representation of an Euclidean vector.
 *
 * This class is a minimum memory container of a vector that offers many operations and tool functions. It is designed
 * to be used as a heap memory structure and is based on std::array for its storage.
 *
 *
 * @tparam Dim The dimension of the vector
 * @tparam TComponent The data type of the vector's components (default to float)
 */
template <size_t Dim, typename TComponent=float>
struct Vector : public std::array<TComponent, Dim>
{
    static constexpr size_t Dimension = Dim;
    using ComponentType = TComponent;

    //////////////////////
    //// Constructors ////
    //////////////////////

    /**
     * Main constructor
     * @param initialize_to_zero If true, initialize the scalar components of the vector to zero
     */
    explicit
    Vector(bool initialize_to_zero = false) {
        if (initialize_to_zero) {
            this->fill(0);
        }
    }

    /**
     * Constructor by a list initializer (ex. Vector<3> A {0, 100, 0}).
     */
    Vector(std::initializer_list<ComponentType> il) {
        std::copy(il.begin(), il.end(), this->begin());
    }

    /**
     * Constructor by c-array or initializer list
     */
    Vector(ComponentType const (&components)[Dimension]) {
        for (size_t i = 0; i < Dimension; ++i) {
            (*this)[i] = components[i];
        }
    }

    /**
     * Constructor by variadic arguments.
     *
     * Example:
     * \code{.cpp}
     * Vector<3> v (0,1,0);
     * \endcode
     *
     * Warning: The data type of the components will the one declared on the class template
     * specifier (float if nothing is specified), hence it will not be the one of the arguments passed.
     *
     * For example, the following vector:
     *
     * \code{.cpp}
     * Vector<3, char> v ((float) 0.2, (float) 1.1, (float) 1.4);
     * \endcode
     *
     * will have its coordinates as char values, and not as float values.
     */
    template<typename... Args>
    Vector(ComponentType first_component, Args&&... other_components)
    {
        recursive_set<0> (first_component, std::forward<Args>(other_components)...);
    }

    ///////////////////
    //// Operators ////
    ///////////////////

    /** Comparison operator with a vector of a same dimension and with component data type of OtherComponentType **/
    template<typename OtherComponentType>
    constexpr bool
    operator==(const Vector<Dimension, OtherComponentType> & other) const
    { return std::equal(this->begin(), this->end(), other.begin()); }

    /** Alias for inner_product **/
    template <size_t OtherDimension, typename OtherComponentType=float>
    constexpr ComponentType
    operator*(const Vector<OtherDimension, OtherComponentType> & other) const
    {
        static_assert(Dimension == OtherDimension, "The inner product cannot be computed on vectors of different size.");
        return inner_product(other);
    }

    /** Alias for cross_product **/
    template <typename OtherComponentType=float>
    constexpr Vector<Dimension, ComponentType>
    operator^(const Vector<Dimension, OtherComponentType> & other) const
    {
        static_assert(Dimension == 3, "The cross product is only implemented for 3-dimensional vectors.");

        return cross_product(other);
    };

    /** Alias for scalar_multiplication **/
    template <class TScalar>
    constexpr Vector<Dimension, ComponentType>
    operator*(const TScalar & scalar) const
    {
        return scalar_multiplication(scalar);
    }

    /** Alias for scalar_division **/
    template <class TScalar>
    constexpr Vector<Dimension, ComponentType>
    operator/(const TScalar & scalar) const
    {
        return scalar_division(scalar);
    }

    /** Alias for direct_sum **/
    template <size_t OtherDimension, typename OtherComponentType=float>
    constexpr Vector<Dimension, ComponentType>
    operator+(const Vector<OtherDimension, OtherComponentType> & other) const
    {
        return direct_sum(other);
    }

    /** Alias for direct_sub **/
    template <size_t OtherDimension, typename OtherComponentType=float>
    constexpr Vector<Dimension, ComponentType>
    operator-(const Vector<OtherDimension, OtherComponentType> & other) const
    {
        return direct_sub(other);
    }

    /////////////////////////////////
    //// Mathematical operations ////
    /////////////////////////////////

    /**
     * Compute the inner product with another vector (aliases are dot, scalar_product and operator*).
     */
    template <typename OtherComponentType=float>
    constexpr ComponentType
    inner_product(const Vector<Dimension, OtherComponentType> & other) const
    {
        return std::inner_product(std::begin(other), std::end(other), std::begin(*this), (TComponent) 0);
    }

    /** Alias for inner_product **/
    template <typename OtherComponentType=float>
    constexpr ComponentType
    dot(const Vector<Dimension, OtherComponentType> & other) const
    {
        return inner_product(other);
    }

    /** Alias for inner_product **/
    template <typename OtherComponentType=float>
    constexpr ComponentType
    scalar_product(const Vector<Dimension, OtherComponentType> & other) const
    {
        return inner_product(other);
    }

    /**
     * Compute the cross product with another vector.
     */
    template <typename OtherComponentType=float>
    constexpr Vector<Dimension, ComponentType>
    cross_product(const Vector<Dimension, OtherComponentType> & other) const
    {
        static_assert(Dimension == 3, "The cross product is only implemented for 3-dimensional vectors.");

        const auto & a = *this;
        const auto & b = other;

        return {
                a[1]*b[2] -a[2]*b[1],
                a[2]*b[0] -a[0]*b[2],
                a[0]*b[1] -a[1]*b[0],
        };
    };

    /** Alias for cross_product **/
    template <typename OtherComponentType=float>
    constexpr Vector<Dimension, ComponentType>
    cross(const Vector<Dimension, OtherComponentType> & other) const
    {
        return cross_product(other);
    };

    /**
     * Compute the direct sum with another vector (sum between each scalar components).
     * @return The resulting vector
     */
    template <typename OtherComponentType=float>
    constexpr Vector<Dimension, ComponentType>
    direct_sum(const Vector<Dimension, OtherComponentType> & other) const
    {
        Vector<Dimension, ComponentType> result(false);
        std::transform(std::begin(other), std::end(other), std::begin(*this), std::begin(result), std::plus<ComponentType >());
        return result;
    }

    /**
     * Compute the direct sub with another vector (subtraction between each scalar components).
     * @return The resulting vector
     */
    template <typename OtherComponentType=float>
    constexpr Vector<Dimension, ComponentType>
    direct_sub(const Vector<Dimension, OtherComponentType> & other) const
    {
        Vector<Dimension, ComponentType> result(false);
        std::transform(std::begin(other), std::end(other), std::begin(*this), std::begin(result), std::minus<ComponentType >());
        return result;
    }

    /**
     * Compute the scalar multiplication (not to be confused by the scalar product).
     * @tparam TScalar The data type of the scalar (ex float, double, int,...)
     * @param scalar The scalar value
     * @return The resulting vector
     */
    template <class TScalar>
    constexpr Vector<Dimension, ComponentType>
    scalar_multiplication(const TScalar & scalar) const
    {
        Vector<Dimension, ComponentType> result(false);
        std::transform(std::begin(*this), std::end(*this), std::begin(result), [scalar] (const ComponentType & component) {
            return component*scalar;
        });
        return result;
    }

    /**
    * Compute the scalar division.
    * @tparam TScalar The data type of the scalar (ex float, double, int,...)
    * @param scalar The scalar value
    * @return The resulting vector
    */
    template <class TScalar>
    constexpr Vector<Dimension, ComponentType>
    scalar_division(const TScalar & scalar) const
    {
        Vector<Dimension, ComponentType> result(false);
        std::transform(std::begin(*this), std::end(*this), std::begin(result), [scalar] (const ComponentType & component) {
            return component/scalar;
        });
        return result;
    }


    /////////////////////////////////
    //// Mathematical properties ////
    /////////////////////////////////

    /** Compute the length of the vector. **/
    constexpr ComponentType
    length() const
    {
        return sqrt((*this).dot(*this));
    };

    /** Get the unit vector (current vector normalized to unit length) **/
    constexpr Vector<Dimension, ComponentType>
    unit() const
    {
        return (*this)/length();
    };


private:
    template<size_t current_index, typename... Args>
    void recursive_set(ComponentType component, Args&&... other_components) {
        (*this)[current_index] = component;
        recursive_set<current_index+1>(std::forward<Args>(other_components)...);
    }

    template<size_t current_index>
    void recursive_set(ComponentType component) {
        (*this)[current_index] = component;
    }
};

template <size_t Dim, typename TComponent=float>
std::ostream& operator<<(std::ostream& os, const Vector<Dim, TComponent>& v)
{
    os << std::string("(");
    os << std::accumulate(std::next(std::begin(v)), std::end(v), std::to_string(v[0]), [](const std::string & s, const TComponent & component) -> std::string {
        return s + std::string(", ") + std::to_string(component);
    });
    os << std::string(")");
    return os;
}

} // namespace algebra

} // namespace caribou

#endif //CARIBOU_ALGEBRA_VECTOR_H
