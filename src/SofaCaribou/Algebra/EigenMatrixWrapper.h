#pragma once
#include <Caribou/macros.h>
#include <Caribou/Traits.h>
#include <sofa/defaulttype/BaseMatrix.h>
#include <Eigen/Core>
#include <Eigen/Sparse>

namespace SofaCaribou::Algebra {

/**
 * The EigenMatrixWrapper is a wrapper around any Eigen matrix. It
 * implements (and overloads) the BaseMatrix class of Sofa. This class hence
 * allows to either map or copy an existng Eigen matrix; and use it in Sofa's
 * components or vistors.
 *
 * \note The wrapper around Eigen sparse matrices has special restrictions.
 *
 * Example:
 * \code{.cpp}
 *    // Wrapper by copy
 *    using EigenDense = Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>;
 *    EigenDense m(100,100);
 *    m(30, 30) = 1;
 *
 *    EigenMatrixWrapper<EigenDense> wrapper(m);
 *    std::cout << m(30, 30) == wrapper(30, 30); // TRUE
 *
 *    wrapper(30, 30) = 100;
 *    std::cout << m(30, 30) == wrapper(30, 30); // FALSE
 * \endcode
 *
 * Example:
 * \code{.cpp}
 *    // Wrapper by reference
 *    using EigenDense = Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>;
 *    EigenDense m(100,100);
 *    m(30, 30) = 1;
 *
 *    EigenMatrixWrapper<EigenDense &> wrapper(m);x
 *    std::cout << m(30, 30) == wrapper(30, 30); // TRUE
 *
 *    wrapper(30, 30) = 100;
 *    std::cout << m(30, 30) == wrapper(30, 30); // TRUE
 * \endcode
 */
template <typename Derived, typename Enable = void>
class EigenMatrixWrapper : public sofa::defaulttype::BaseMatrix
{
    static_assert(
        std::is_base_of_v<Eigen::EigenBase<std::decay_t<Derived> >, std::decay_t<Derived> >,
        "The class template argument 'Derived' must be derived from Eigen::EigenBase<Derived>."
    );

public:
    using EigenType = std::remove_cv_t<std::remove_reference_t<Derived>>;
    using Base = sofa::defaulttype::BaseMatrix;
    using Index = Base::Index;
    using Real = SReal;

    EigenMatrixWrapper(std::remove_reference_t<Derived> & eigen_matrix) : p_eigen_matrix(eigen_matrix) {}

    // Abstract methods overrides
    inline Index rowSize() const final { return p_eigen_matrix.rows(); }
    inline Index colSize() const final { return p_eigen_matrix.cols(); }
    inline Real  element(Index i, Index j) const final { return p_eigen_matrix(i,j); }

    /** Resize the matrix to nbRow x nbCol dimensions. This method resets to zero all entries. */
    inline void  resize(Index nbRow, Index nbCol) final { p_eigen_matrix.setConstant(nbRow, nbCol, 0); }

    /** Set all entries to zero. Keeps the current matrix dimensions. */
    inline void  clear() override {p_eigen_matrix.setZero();}

    /** Set this value of the matrix entry (i, j) to the value of v. */
    inline void  set(Index i, Index j, double v) final {
        p_eigen_matrix(i, j) = v;
    }

    /** Adds v to the value of the matrix entry (i, j). */
    inline void  add(Index i, Index j, double v) final {
        p_eigen_matrix(i, j) += v;
    }

    // Block operations on 3x3 and 2x2 sub-matrices
    inline void add(Index i, Index j, const sofa::defaulttype::Mat3x3d & m) override {return add_block(i, j, m);}
    inline void add(Index i, Index j, const sofa::defaulttype::Mat3x3f & m) override {return add_block(i, j, m);}
    inline void add(Index i, Index j, const sofa::defaulttype::Mat2x2d & m) override {return add_block(i, j, m);}
    inline void add(Index i, Index j, const sofa::defaulttype::Mat2x2f & m) override {return add_block(i, j, m);}

    /** Sets the entire row i to zero */
    inline void clearRow(Index i) final {
        p_eigen_matrix.row(i) *= 0;
    }

    /** Sets the rows from index imin to index imax (inclusively) to zero */
    inline void clearRows(Index imin, Index imax) final {
        p_eigen_matrix.middleRows(imin, (imax-imin)+1) *= 0;
    }

    /** Sets the entire columns i to zero */
    inline void clearCol(Index i) final {
        p_eigen_matrix.col(i) *= 0;
    }

    /** Sets the columns from index imin to index imax (inclusively) to zero */
    inline void clearCols(Index imin, Index imax) final {
        p_eigen_matrix.middleCols(imin, (imax-imin)+1) *= 0;
    }

private:

    template <typename Scalar, std::size_t nrows, std::size_t ncols>
    void add_block(Index i, Index j, const sofa::defaulttype::Mat<nrows, ncols, Scalar> & m) {
        const Eigen::Map<const Eigen::Matrix<Scalar, nrows, ncols, Eigen::RowMajor>> block(&(m[0][0]));
        p_eigen_matrix.block(i, j, nrows, ncols) += block. template cast<typename EigenType::Scalar>();
    }

    ///< The actual Eigen Matrix.
    Derived p_eigen_matrix;
};

///////////////////////////////////
/// SparseMatrix specialization ///
///////////////////////////////////
template <typename Derived>
class EigenMatrixWrapper<Derived, CLASS_REQUIRES(std::is_base_of_v<Eigen::SparseMatrixBase<std::decay_t<Derived>>, std::decay_t<Derived>>)> : public sofa::defaulttype::BaseMatrix
{

public:
    using EigenType = std::remove_cv_t<std::remove_reference_t<Derived>>;
    using Base = sofa::defaulttype::BaseMatrix;
    using Index = Base::Index;
    using Real = SReal;

    EigenMatrixWrapper(std::remove_reference_t<Derived> & eigen_matrix) : p_eigen_matrix(eigen_matrix) {}

    // Abstract methods overrides
    inline Index rowSize() const final { return p_eigen_matrix.rows(); }
    inline Index colSize() const final { return p_eigen_matrix.cols(); }
    /**
     * @brief Return the matrix entry (i,j).
     * \warning If the matrix hasn't been initialized by calling compress() or set(), this
     *          method will always return 0.
     */
    inline Real  element(Index i, Index j) const final {
        caribou_assert(p_initialized && "Accessing an element on an uninitialized matrix.");
        return this->p_eigen_matrix.coeff(i,j);
    }

    /** Resize the matrix to nbRow x nbCol dimensions. This method resets to zero all entries. */
    inline void  resize(Index nbRow, Index nbCol) final {
        p_triplets.clear();
        this->p_eigen_matrix.resize(nbRow, nbCol);
        p_initialized = false;
    }

    /** Set all entries to zero. Keeps the current matrix dimensions. */
    inline void  clear() final {
        p_triplets.clear();
        this->p_eigen_matrix.setZero();
        p_initialized = false;
    }

    /**
     * Set this value of the matrix entry (i, j) to the value of v.
     *
     * \warning When the matrix hasn't been initialized yet, this method will have to do it. This
     *          means that further calls to set or add will be much slower and generate an uncompressed
     *          matrix. If you have to use the set method, make sure that all calls to add has been
     *          made beforehand.
     */
    inline void  set(Index i, Index j, double v) final {
        if (not p_initialized) {
            initialize();
        }

        this->p_eigen_matrix.coeffRef(i, j) = v;
    }

    /**
     * Adds v to the value of the matrix entry (i, j).
     *
     * \warning If this method is called after the matrix has been initialized, the matrix will
     *          become uncompressed.
     */
    inline void  add(Index i, Index j, double v) final {
        // Note: this "if" condition should't slow down that much the addition of multiple entries
        //       because of branch predictions (the conditional branch will be same for the next
        //       X calls to add until compress is called).
        if (not p_initialized) {
            p_triplets.emplace_back(i, j, v);
        } else {
            p_eigen_matrix.coeffRef(i, j) += v;
        }
    }

    /**
     * Compress the matrix.
     *
     * If it is the first time that this method is called, then the matrix is built from a list of
     * triplets (i, j, value) accumulated during the calls of the methodes 'add' and 'set'.
     */
    void compress() final {
        if (not p_initialized) {
            initialize();
        } else {
            p_eigen_matrix.makeCompressed();
        }
    }

    // Block operations on 3x3 and 2x2 sub-matrices
    inline void add(Index i, Index j, const sofa::defaulttype::Mat3x3d & m) override {return add_block(i, j, m);}
    inline void add(Index i, Index j, const sofa::defaulttype::Mat3x3f & m) override {return add_block(i, j, m);}
    inline void add(Index i, Index j, const sofa::defaulttype::Mat2x2d & m) override {return add_block(i, j, m);}
    inline void add(Index i, Index j, const sofa::defaulttype::Mat2x2f & m) override {return add_block(i, j, m);}

    /** Sets the entire row i to zero */
    inline void clearRow(Index i) final {
        if (not p_initialized) {
            initialize();
        }

#if EIGEN_VERSION_AT_LEAST(3,3,0)
        p_eigen_matrix.row(i) *= static_cast<const typename EigenType::Scalar &>(0);
#else
        p_eigen_matrix.middleRows(i,1) *= static_cast<const typename EigenType::Scalar &>(0);
#endif
    }

    /** Sets the rows from index imin to index imax (inclusively) to zero */
    inline void clearRows(Index imin, Index imax) final {
        if (not p_initialized) {
            initialize();
        }

        p_eigen_matrix.middleRows(imin, (imax-imin)+1) *= static_cast<const typename EigenType::Scalar &>(0);
    }

    /** Sets the entire columns i to zero */
    inline void clearCol(Index i) final {
        if (not p_initialized) {
            initialize();
        }

        p_eigen_matrix.col(i) *= static_cast<const typename EigenType::Scalar &>(0);
    }

    /** Sets the columns from index imin to index imax (inclusively) to zero */
    inline void clearCols(Index imin, Index imax) final {
        if (not p_initialized) {
            initialize();
        }

        p_eigen_matrix.middleCols(imin, (imax-imin)+1) *= static_cast<const typename EigenType::Scalar &>(0);
    }

private:

    /**
     * Block addition with a NxC matrix
     */
    template <typename Scalar, std::size_t N, std::size_t C>
    void add_block(Index i, Index j, const sofa::defaulttype::Mat<N, C, Scalar> & m) {
        for (size_t k=0;k<N;++k) {
            for (size_t l=0;l<C;++l) {
                const auto value = static_cast<typename EigenType::Scalar>(m[k][l]);
                if (not p_initialized) {
                    p_triplets.emplace_back(i+k, j+l, value);
                } else {
                    p_eigen_matrix.coeffRef(i+k, j+l) += value;
                }
            }
        }
    }

    /**
     * @brief initialize the matrix with the accumulated triplets. The resulting matrix will be
     *        compressed. Once this is done, the operation add and addblock will be much slower
     *        if the entries to be added were zero before hand.
     */
    void initialize() {
        p_eigen_matrix.setFromTriplets(p_triplets.begin(), p_triplets.end());
        p_triplets.clear();
        p_initialized = true;
    }

    ///< Triplets are used to store matrix entries before the call to 'compress'.
    /// Duplicates entries are summed up.
    std::vector<Eigen::Triplet<typename EigenType::Scalar>> p_triplets;

    ///< Whether or not the matrix has been initialized with triplets yet. This will determined the behavior
    /// of the add and set methods. When the matrix hasn't been initialized, the add and set methods will simply append
    /// the new coefficient to the triplet list.
    bool p_initialized = false;

    ///< The actual Eigen Matrix.
    Derived p_eigen_matrix;
};

} // namespace SofaCaribou::Algebra
