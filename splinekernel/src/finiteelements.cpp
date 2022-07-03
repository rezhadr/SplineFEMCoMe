#include "finiteelements.hpp"
#include "basisfunctions.hpp"
#include "utilities.hpp"

#include <cmath>
#include <algorithm>
#include <numeric>

namespace cie
{
namespace splinekernel
{

// This function evaluates the one-dimensional B-Spline basis functions 
// or their first derivative. Should have been in basisfunctions.cpp,
// but we forgot :)
double evaluateBSplineDerivative( double t, 
                                  size_t functionIndex, 
                                  size_t degree, 
                                  const std::vector<double>& knotVector, 
                                  size_t diffOrder )
{
    if( diffOrder == 0 )
    {
        return evaluateBSplineBasis( t, functionIndex, degree, knotVector );
    }
    else if( diffOrder == 1 )
    {
        if( degree == 0 )
        {
            throw std::runtime_error( "Invalid polynomial degree!" );
        }

        double factor1 = knotVector[functionIndex + degree] - knotVector[functionIndex];
        double factor2 = knotVector[functionIndex + degree + 1] - knotVector[functionIndex + 1];

        if( std::abs( factor1 ) > 1e-10 )
        {
            factor1 = degree / factor1 * evaluateBSplineBasis( t, functionIndex, degree - 1, knotVector );
        }

        if( std::abs( factor2 ) > 1e-10 )
        {
            factor2 = degree / factor2 * evaluateBSplineBasis( t, functionIndex + 1, degree - 1, knotVector );
        }

        return factor1 - factor2;
    }
    else
    {
        throw std::runtime_error( "Invalid diff order." );
    }
}

namespace detail
{

std::array<double, 2> mapToGlobalCoordinates( std::array<double, 2> localCoordinates,
                                              std::array<size_t, 2> elementIndices,
                                              std::array<double, 2> lengths,
                                              std::array<double, 2> origin,
                                              std::array<size_t, 2> numberOfElements )
{
    runtime_check(numberOfElements[0] != 0 && numberOfElements[1] != 0,
                  "Number of elements cannot be zero");

    std::array<double, 2> globCoor;

    for (size_t i = 0; i < 2; i++)
    {
        globCoor[i] = (elementIndices[i] + (localCoordinates[i] + 1)/ 2) *
                       lengths[i] / numberOfElements[i] +
                       origin[i];
    }
    
    return globCoor;
}

KnotVectors constructOpenKnotVectors( std::array<size_t, 2> numberOfElements,
                                      std::array<size_t, 2> polynomialDegrees,
                                      std::array<size_t, 2> continuities,
                                      std::array<double, 2> lengths,
                                      std::array<double, 2> origin )
{
    runtime_check(polynomialDegrees[0] > continuities[0] && polynomialDegrees[1] > continuities[1],
                  "Continuity must be lower than polynomial degree");

    KnotVectors openKnotVectors = {std::vector<double>(polynomialDegrees[0]+1,origin[0]),
                                   std::vector<double>(polynomialDegrees[1]+1,origin[1])};

    for (size_t i = 0; i < 2; i++)
    {
        for (size_t j = 1; j < numberOfElements[i]; j++)
        {
            for (size_t k = 0; k < polynomialDegrees[i]-continuities[i]; k++)
            {
                openKnotVectors[i].push_back(
                    double(j) / numberOfElements[i] * lengths[i] + origin[i]
                );
            }
        }

        for (size_t k = 0; k < polynomialDegrees[i]+1; k++)
            {
                openKnotVectors[i].push_back(
                    lengths[i] + origin[i]
                );
            }   
    }
    
    return openKnotVectors;
}

size_t findKnotSpan( double min, double max, size_t n, double x )
{
    runtime_check(n > 0, "Number of span is at least 1");

    if (std::abs(x - min) < 1e-12)
        return 0;

    if (std::abs(x - max) < 1e-12)
        return n - 1;
    
    return std::floor((x-min)/(max-min)*n);
}

} // namespace detail

BSplineFiniteElementPatch::BSplineFiniteElementPatch( std::array<size_t, 2> numberOfElements,
                                                      std::array<size_t, 2> polynomialDegrees,
                                                      std::array<size_t, 2> continuities,
                                                      std::array<double, 2> lengths,
                                                      std::array<double, 2> origin,
                                                      IntegrationPointProvider integrationPointProvider ) :
                                                      numberOfElements_(numberOfElements),
                                                      polynomialDegrees_(polynomialDegrees),
                                                      continuities_(continuities),
                                                      lengths_(lengths),
                                                      origin_(origin),
                                                      integrationPointProvider_(integrationPointProvider)
{
    knotVectors_ = detail::constructOpenKnotVectors(
        numberOfElements_,
        polynomialDegrees_,
        continuities_,
        lengths_,
        origin_
    );
} // constructor

// When we call this function internally we often know the kont span/element indices in which
// globalCoordinates lie already. But this being a public function it would make little sense to
// ask a user to compute the knot span indices which he shouldn't need to care about. 
std::vector<double> BSplineFiniteElementPatch::evaluateActiveBasisAt( std::array<double, 2> globalCoordinates,
                                                                      std::array<size_t, 2> diffOrders ) const
{
    std::array<size_t, 2> knotSpan = {
        detail::findKnotSpan(
        origin_[0], origin_[0]+lengths_[0], numberOfElements_[0], globalCoordinates[0]
        ),
        detail::findKnotSpan(
        origin_[1], origin_[1]+lengths_[1], numberOfElements_[1], globalCoordinates[1]
        )};
    
    std::array<size_t, 2> shapeInd = {
        knotSpan[0] * (polynomialDegrees_[0] - continuities_[0]),
        knotSpan[1] * (polynomialDegrees_[1] - continuities_[1])
    };

    std::vector<double> basisFunc;

    for (size_t i = 0; i < polynomialDegrees_[0]+1; i++)
    {
        for (size_t j = 0; j < polynomialDegrees_[1]+1; j++)
        {
            basisFunc.push_back(
                evaluateBSplineDerivative(
                    globalCoordinates[0],
                    shapeInd[0] + i,
                    polynomialDegrees_[0],
                    knotVectors_[0],
                    diffOrders[0]
                ) *
                evaluateBSplineDerivative(
                    globalCoordinates[1],
                    shapeInd[1] + j,
                    polynomialDegrees_[1],
                    knotVectors_[1],
                    diffOrders[1]
                )
            );
        }
    }
    
    return basisFunc;
}

ElementLinearSystem BSplineFiniteElementPatch::integrateElementSystem( std::array<size_t, 2> elementIndices,
                                                                       const SpatialFunction& sourceFunction ) const
{
    size_t nDOFs = (polynomialDegrees_[0] + 1) * (polynomialDegrees_[1] + 1);

    linalg::Matrix K(nDOFs, nDOFs, 0.0);
    std::vector<double> F(nDOFs, 0.0);

    IntegrationPoints iPointXi = integrationPointProvider_(polynomialDegrees_[0]+1);
    IntegrationPoints iPointEt = integrationPointProvider_(polynomialDegrees_[1]+1);

    for (size_t i = 0; i < iPointXi[0].size(); i++)
    {
        for (size_t j = 0; j < iPointEt[0].size(); j++)
        {
            std::array<double, 2> locCoor = {iPointXi[0].at(i), iPointEt[0].at(j)};

            double Wij = iPointXi[1].at(i) * iPointEt[1].at(j);

            auto globCoor = detail::mapToGlobalCoordinates(
                locCoor, elementIndices, lengths_, origin_, numberOfElements_
            );

            auto N = evaluateActiveBasisAt(
                globCoor, {0,0}
            );

            auto dNdx = evaluateActiveBasisAt(
                globCoor, {1,0}
            );

            auto dNdy = evaluateActiveBasisAt(
                globCoor, {0,1}
            );

            double detJ = lengths_[0] / numberOfElements_[0] *
                          lengths_[1] / numberOfElements_[1] / 4;

            auto f = sourceFunction(globCoor[0], globCoor[1]);

            for (size_t p = 0; p < nDOFs; p++)
            {
                for (size_t q = 0; q < nDOFs; q++)
                {
                    K(p,q) +=  (
                        dNdx[p] * dNdx[q] +
                        dNdy[p] * dNdy[q]
                    ) * Wij * detJ;
                }
                
                F[p] += N[p] * f * Wij * detJ;
            }
            
        }
        
    }
    
    return {K, F};
}

} // namespace splinekernel
} // namespace cie
