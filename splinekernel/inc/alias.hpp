#pragma once

#include <algorithm>

namespace cie
{
namespace splinekernel
{
    class CompressedSparseRowMatrix;

    using GlobalLinearSystem = std::pair<CompressedSparseRowMatrix, std::vector<double>>;
    using ElementLinearSystem = std::pair<linalg::Matrix, std::vector<double>>;

    using SpatialFunction = std::function<double( double, double )>;

    using IntegrationPoints = std::array<std::vector<double>, 2>;
    using IntegrationPointProvider = std::function<IntegrationPoints( size_t )>;

    using KnotVectors = std::array<std::vector<double>, 2>;
    using LocationMap = std::vector<size_t>;
}
}
