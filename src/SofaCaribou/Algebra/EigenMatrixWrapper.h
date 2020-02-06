#pragma once
#include <Caribou/Traits.h>
#include <sofa/defaulttype/BaseMatrix.h>
#include <Eigen/Core>
#include <Eigen/Sparse>

namespace SofaCaribou::Algebra {

/**
 * The EigenMatrixWrapper is a wrapper around any Eigen matrix types that
 * implements (and overloads) the BaseMatrix class of Sofa.
 *
 * @example
 * \code{.cpp}
 *
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

    template <typename Scalar, int nrows, int ncols>
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
     * \warning When the matrix hasn't been compressed yet, this method will have to do it. This
     *          means that further calls to set or add will be much slower. If you have to use the
     *          set method, make sure that all calls to add has been made beforehand.
     */
    inline void  set(Index i, Index j, double v) final {
        if (not p_initialized) {
            p_eigen_matrix.setFromTriplets(p_triplets.begin(), p_triplets.end());
            p_triplets.clear();
            p_initialized = true;
        }

        this->p_eigen_matrix.coeffRef(i, j) = v;
    }

    /** Adds v to the value of the matrix entry (i, j). */
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
            p_eigen_matrix.setFromTriplets(p_triplets.begin(), p_triplets.end());
            p_triplets.clear();
            p_initialized = true;
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

    /**
     * Block addition with a NxC matrix
     */
    template <typename Scalar, int N, int C>
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
