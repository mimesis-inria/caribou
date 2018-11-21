#ifndef CARIBOU_ALGEBRA_VECTOR_H
#define CARIBOU_ALGEBRA_VECTOR_H

#include <ostream>
#include <cstddef>
#include <array>
#include <initializer_list>
#include <algorithm>
#include <numeric>
#include <cmath>

#include <Caribou/config.h>

namespace caribou
{
namespace algebra
{

/**
 * A simple representation of a vector.
 *
 * This class is a minimum memory container of a vector that offers many operations and tool functions. It is designed
 * to be used as a heap memory structure and is based on std::array for its storage.
 *
 *
 * @tparam Dim The dimension of the vector
 * @tparam TValueType The data type of the vector's components (default to float)
 */
template <size_t Dim, typename TValueType=FLOATING_POINT_TYPE>
struct Vector : public std::array<TValueType, Dim>
{
    static constexpr size_t Dimension = Dim;
    using ValueType = TValueType;

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
    template <typename OtherValueType>
    Vector(std::initializer_list<OtherValueType> il) {
        std::copy(il.begin(), il.end(), this->begin());
    }

    /**
     * Constructor by c-array or initializer list
     */
    template <typename OtherValueType>
    Vector(OtherValueType const (&components)[Dimension]) {
        for (size_t i = 0; i < Dimension; ++i) {
            (*this)[i] = static_cast<ValueType>(components[i]);
        }
    }

    /**
     * Copy constructor from another vector of a different data type
     */
    template <typename OtherValueType>
    Vector(const Vector<Dimension, OtherValueType> & other)
    {
         for (size_t i = 0; i < Dimension; ++i) {
             (*this)[i] = static_cast<ValueType> (other[i]);
         }
    }

    /**
     * Copy constructor from another type of vector (it must implement the [] operator and its data type must match ValueType)
     */
    template <typename OtherVectorType>
    Vector(const OtherVectorType & other)
    {
        for (size_t i = 0; i < Dimension; ++i) {
            (*this)[i] = static_cast<ValueType> (other[i]);
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
    Vector(ValueType first_component, Args&&... other_components)
    {
        recursive_set<0> (first_component, std::forward<Args>(other_components)...);
    }

    ///////////////////
    //// Operators ////
    ///////////////////

    /** Comparison operator with a vector of a same dimension and with component data type of OtherValueType **/
    template<typename OtherValueType>
    constexpr bool
    operator==(const Vector<Dimension, OtherValueType> & other) const
    { return std::equal(this->begin(), this->end(), other.begin()); }

    /** Alias for inner_product **/
    template <typename OtherValueType>
    constexpr ValueType
    operator*(const Vector<Dimension, OtherValueType> & other) const
    {
        return inner_product(other);
    }

    /** Alias for cross_product **/
    template <typename OtherValueType>
    constexpr Vector<Dimension, ValueType>
    operator^(const Vector<Dimension, OtherValueType> & other) const
    {
        static_assert(Dimension == 3, "The cross product is only implemented for 3-dimensional vectors.");

        return cross_product(other);
    }

    /** Alias for scalar_multiplication **/
    template <class TScalar>
    constexpr Vector<Dimension, ValueType>
    operator*(const TScalar & scalar) const
    {
        return scalar_multiplication(scalar);
    }

    /** Alias for scalar_division **/
    template <class TScalar>
    constexpr Vector<Dimension, ValueType>
    operator/(const TScalar & scalar) const
    {
        return scalar_division(scalar);
    }

    /** Alias for direct_sum **/
    template <typename OtherValueType>
    constexpr Vector<Dimension, ValueType>
    operator+(const Vector<Dimension, OtherValueType> & other) const
    {
        return direct_sum(other);
    }

    /** Alias for direct_sub **/
    template <typename OtherValueType=float>
    constexpr Vector<Dimension, ValueType>
    operator-(const Vector<Dimension, OtherValueType> & other) const
    {
        return direct_sub(other);
    }

    /////////////////////////////////
    //// Mathematical operations ////
    /////////////////////////////////

    /**
     * Compute the inner product with another vector (aliases are dot, scalar_product and operator*).
     */
    template <typename OtherValueType>
    constexpr ValueType
    inner_product(const Vector<Dimension, OtherValueType> & other) const
    {
        return std::inner_product(std::begin(other), std::end(other), std::begin(*this), (ValueType) 0);
    }

    /** Alias for inner_product **/
    template <typename OtherValueType>
    constexpr ValueType
    dot(const Vector<Dimension, OtherValueType> & other) const
    {
        return inner_product(other);
    }

    /** Alias for inner_product **/
    template <typename OtherValueType>
    constexpr ValueType
    scalar_product(const Vector<Dimension, OtherValueType> & other) const
    {
        return inner_product(other);
    }

    /**
     * Compute the cross product with another vector.
     */
    template <typename OtherValueType>
    constexpr Vector<Dimension, ValueType>
    cross_product(const Vector<Dimension, OtherValueType> & other) const
    {
        static_assert(Dimension == 3, "The cross product is only implemented for 3-dimensional vectors.");

        const auto & a = *this;
        const auto & b = other;

        return {
                a[1]*b[2] -a[2]*b[1],
                a[2]*b[0] -a[0]*b[2],
                a[0]*b[1] -a[1]*b[0],
        };
    }

    /** Alias for cross_product **/
    template <typename OtherValueType>
    constexpr Vector<Dimension, ValueType>
    cross(const Vector<Dimension, OtherValueType> & other) const
    {
        return cross_product(other);
    }

    /**
     * Compute the direct sum with another vector (sum between each scalar components).
     * @return The resulting vector
     */
    template <typename OtherValueType>
    constexpr Vector<Dimension, ValueType>
    direct_sum(const Vector<Dimension, OtherValueType> & other) const
    {
        Vector<Dimension, ValueType> result(false);
        std::transform(std::begin(other), std::end(other), std::begin(*this), std::begin(result), std::plus<ValueType >());
        return result;
    }

    /**
     * Compute the direct sub with another vector (subtraction between each scalar components).
     * @return The resulting vector
     */
    template <typename OtherValueType>
    constexpr Vector<Dimension, ValueType>
    direct_sub(const Vector<Dimension, OtherValueType> & other) const
    {
        Vector<Dimension, ValueType> result(false);
        std::transform(std::begin(other), std::end(other), std::begin(*this), std::begin(result), std::minus<ValueType >());
        return result;
    }

    /**
     * Compute the direct multiplication with another vector (multiplication between each scalar components).
     * @return The resulting vector
     */
    template <typename OtherValueType>
    constexpr Vector<Dimension, ValueType>
    direct_mult(const Vector<Dimension, OtherValueType> & other) const
    {
        Vector<Dimension, ValueType> result(false);
        std::transform(std::begin(other), std::end(other), std::begin(*this), std::begin(result), std::multiplies<ValueType >());
        return result;
    }

    /**
     * Compute the direct division with another vector (division between each scalar components).
     * @return The resulting vector
     */
    template <typename OtherValueType>
    constexpr Vector<Dimension, ValueType>
    direct_division(const Vector<Dimension, OtherValueType> & other) const
    {
        Vector<Dimension, ValueType> result(false);
        std::transform(std::begin(*this), std::end(*this), std::begin(other), std::begin(result), std::divides<ValueType >());
        return result;
    }

    /**
     * Compute the scalar multiplication (not to be confused by the scalar product).
     * @tparam TScalar The data type of the scalar (ex float, double, int,...)
     * @param scalar The scalar value
     * @return The resulting vector
     */
    template <class TScalar>
    constexpr Vector<Dimension, ValueType>
    scalar_multiplication(const TScalar & scalar) const
    {
        Vector<Dimension, ValueType> result(false);
        std::transform(std::begin(*this), std::end(*this), std::begin(result), [scalar] (const ValueType & component) {
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
    constexpr Vector<Dimension, ValueType>
    scalar_division(const TScalar & scalar) const
    {
        Vector<Dimension, ValueType> result(false);
        std::transform(std::begin(*this), std::end(*this), std::begin(result), [scalar] (const ValueType & component) {
            return component/scalar;
        });
        return result;
    }


    /////////////////////////////////
    //// Mathematical properties ////
    /////////////////////////////////

    /** Compute the length of the vector. **/
    constexpr ValueType
    length() const
    {
        return sqrt((*this).dot(*this));
    };

    /** Compute the length^2 of the vector. **/
    constexpr ValueType
    length_squared() const
    {
        return (*this).dot(*this);
    }

    /** Get the unit vector (current vector normalized to unit length) **/
    constexpr Vector<Dimension, ValueType>
    unit() const
    {
        return (*this)/length();
    }

    /////////////////////////////////
    ////      Tool functions     ////
    /////////////////////////////////
    inline std::string
    to_string() const
    {
        return std::string("(") +
               std::accumulate(std::next(this->begin()), this->end(), std::to_string(this->at(0)), [](const std::string & s, const ValueType & component) -> std::string {
                   return s + std::string(", ") + std::to_string(component);
               }) +
               std::string(")");
    }


private:
    template<size_t current_index, typename... Args>
    void recursive_set(ValueType component, Args&&... other_components) {
        (*this)[current_index] = component;
        recursive_set<current_index+1>(std::forward<Args>(other_components)...);
    }

    template<size_t current_index>
    void recursive_set(ValueType component) {
        (*this)[current_index] = component;
    }
};

} // namespace algebra

} // namespace caribou

template <size_t Dim, typename TComponent=FLOATING_POINT_TYPE>
inline std::ostream&
operator<<(std::ostream& os, const caribou::algebra::Vector<Dim, TComponent>& v)
{
    os << std::string("(");
    os << std::accumulate(std::next(std::begin(v)), std::end(v), std::to_string(v[0]), [](const std::string & s, const TComponent & component) -> std::string {
        return s + std::string(", ") + std::to_string(component);
    });
    os << std::string(")");
    return os;
}

template <size_t Dim, typename TOtherType, typename TComponent=FLOATING_POINT_TYPE>
inline caribou::algebra::Vector<Dim, TComponent>
operator * (const TOtherType & a, const caribou::algebra::Vector<Dim, TComponent> & v)
{
    return v*a;
}

#endif //CARIBOU_ALGEBRA_VECTOR_H
