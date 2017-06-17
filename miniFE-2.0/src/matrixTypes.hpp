#ifndef minife_matrix_types_hpp

#include <elementTypes.hpp>
#ifdef MINIFE_ELL_MATRIX
#include <ELLMatrix.hpp>
typedef miniFE::ELLMatrix MatrixType;
#else
#include <CSRMatrix.hpp>
typedef miniFE::CSRMatrix MatrixType;
#endif
#include <Vector.hpp>
typedef miniFE::Vector VectorType;

#endif
