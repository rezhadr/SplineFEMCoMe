#include "linalg.hpp"

namespace cie
{
namespace splinekernel
{

using VectorOfMatrices = std::vector<linalg::Matrix>;

VectorOfMatrices evaluateSurface( const std::array<std::vector<double>, 2>& knotVectors,
                                  const VectorOfMatrices& controlPoints,
                                  std::array<size_t, 2> numberOfSamplePoints );

} // splinekernel
} // cie
