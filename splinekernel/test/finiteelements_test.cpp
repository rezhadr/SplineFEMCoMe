#include "catch.hpp"
#include "finiteelements.hpp"

#include <vector>

namespace cie
{
namespace splinekernel
{


TEST_CASE( "BSplineFiniteElementPatch_mapToGlobalCoordinates_test" )
{
    std::array<double, 2> lengths { 3.0, 4.5 };
    std::array<double, 2> origin { -1.5, 2.5 };
    std::array<size_t, 2> numberOfElements { 2, 3 };

    std::array<double, 2> rs1, rs2, rs3, rs4, rs5, rs6;

    // Valid cases
    REQUIRE_NOTHROW( rs1 = detail::mapToGlobalCoordinates( {  0.8,  0.2 }, { 0, 0 }, lengths, origin, numberOfElements ) );
    REQUIRE_NOTHROW( rs2 = detail::mapToGlobalCoordinates( { -0.3,  0.4 }, { 0, 1 }, lengths, origin, numberOfElements ) );
    REQUIRE_NOTHROW( rs3 = detail::mapToGlobalCoordinates( { -1.0,  0.6 }, { 0, 2 }, lengths, origin, numberOfElements ) );
    REQUIRE_NOTHROW( rs4 = detail::mapToGlobalCoordinates( { -0.9, -0.1 }, { 1, 0 }, lengths, origin, numberOfElements ) );
    REQUIRE_NOTHROW( rs5 = detail::mapToGlobalCoordinates( {  0.7, -0.5 }, { 1, 1 }, lengths, origin, numberOfElements ) );
    REQUIRE_NOTHROW( rs6 = detail::mapToGlobalCoordinates( {  0.3,  0.9 }, { 1, 2 }, lengths, origin, numberOfElements ) );

    CHECK( rs1[0] == Approx( -0.15 ) );
    CHECK( rs1[1] == Approx( 3.4 ) );

    CHECK( rs2[0] == Approx( -0.975 ) );
    CHECK( rs2[1] == Approx( 5.05 ) );

    CHECK( rs3[0] == Approx( -1.5 ) );
    CHECK( rs3[1] == Approx( 6.7 ) );
        
    CHECK( rs4[0] == Approx( 0.075 ) );
    CHECK( rs4[1] == Approx( 3.175 ) );

    CHECK( rs5[0] == Approx( 1.275 ) );
    CHECK( rs5[1] == Approx( 4.375 ) );
        
    CHECK( rs6[0] == Approx( 0.975 ) );
    CHECK( rs6[1] == Approx( 6.925 ) );

    // Degenerate cases with at least one number of elements being zero
    CHECK_THROWS( detail::mapToGlobalCoordinates( { 0.5, 0.5 }, { 0, 0 }, lengths, origin, { 0, 0 } ) );
    CHECK_THROWS( detail::mapToGlobalCoordinates( { 0.5, 0.5 }, { 0, 0 }, lengths, origin, { 0, 1 } ) );
    CHECK_THROWS( detail::mapToGlobalCoordinates( { 0.5, 0.5 }, { 0, 0 }, lengths, origin, { 1, 0 } ) );

} // BSplineFiniteElementPatch_mapToGlobalCoordinates_test



TEST_CASE( "BSplineFiniteElementPatch_constructOpenKnotVector_test" )
{
    // Valid case
    auto result = detail::constructOpenKnotVectors( { 5, 7 }, { 2, 3 }, { 1, 2 }, { 3, 5 }, { -2, 4 } );

    KnotVectors expected = 
    { {
        { -2.0, -2.0, -2.0, -1.4, -0.8, -0.2, 0.4, 1.0, 1.0, 1.0 },
        {  4.0, 4.0, 4.0, 4.0, 4.71429, 5.42857, 6.14286, 6.85714, 7.57143, 8.28571, 9.0, 9.0, 9.0, 9.0 } 
    } };

    for( size_t axis = 0; axis < expected.size( ); ++axis )
    {
        REQUIRE( result[axis].size( ) == expected[axis].size( ) );

        for( size_t index = 0; index < expected[axis].size( ); ++index )
        {
            CHECK( result[axis][index] == Approx( expected[axis][index] ) );
        }
    }

    // Degenerate case: Continuity is higher than the polynomial degree
    std::array<size_t, 2> polynomialDegrees = { 2, 2 };

    std::array<size_t, 2> continuities1 = { 2, 1 };
    std::array<size_t, 2> continuities2 = { 1, 2 };
    std::array<size_t, 2> continuities3 = { 2, 2 };

    CHECK_THROWS( detail::constructOpenKnotVectors( { 2, 3 }, polynomialDegrees, continuities1, { 1, 1 }, { 1, 1 } ) );
    CHECK_THROWS( detail::constructOpenKnotVectors( { 2, 3 }, polynomialDegrees, continuities2, { 1, 1 }, { 1, 1 } ) );
    CHECK_THROWS( detail::constructOpenKnotVectors( { 2, 3 }, polynomialDegrees, continuities3, { 1, 1 }, { 1, 1 } ) );
}



TEST_CASE( "findKnotSpan_test" )
{
    // Spans: [ -2.0, -1.143, -0.286,  0.571,  1.429, 2.286, 3.143, 4.0 ]
    CHECK( detail::findKnotSpan( -2.0, 4.0, 7, -2.0 ) == 0 );
    CHECK( detail::findKnotSpan( -2.0, 4.0, 7, 0.570 ) == 2 );
    CHECK( detail::findKnotSpan( -2.0, 4.0, 7, 0.572 ) == 3 );
    CHECK( detail::findKnotSpan( -2.0, 4.0, 7, 4.0 ) == 6 );
	
	CHECK_THROWS( detail::findKnotSpan( -2.0, 4.0, 0, -2.0 ) );
}



TEST_CASE( "BSplineFiniteElementPatch_constructor_test" )
{
    auto dummy_integrator = []( size_t )
    {
        return IntegrationPoints { };
    };

    CHECK_NOTHROW( BSplineFiniteElementPatch( { 2, 3 }, { 3, 2 }, { 2, 1 }, { 1.0, 1.0 }, { 0.0, 0.0 }, dummy_integrator ) );
    CHECK_NOTHROW( BSplineFiniteElementPatch( { 2, 3 }, { 1, 1 }, { 0, 0 }, { 1.0, 1.0 }, { 1.0, 1.0 }, dummy_integrator ) );
    CHECK_NOTHROW( BSplineFiniteElementPatch( { 1, 1 }, { 2, 3 }, { 1, 1 }, { 1.0, 1.0 }, { 1.0, 1.0 }, dummy_integrator ) );
}



TEST_CASE( "BSplineFiniteElementPatch_evaluateActiveBasisAt_test" )
{
    // Expected tensor product values for three evaluation points and no derivatives: Ni(r) * Nj(s)
    std::vector<std::vector<double>> expectedBasis_00 = 
    {
        { 1.20000000e-02, 1.15703704e-01, 9.21481481e-02, 2.37037037e-03, 3.90000000e-02, 3.76037037e-01,
          2.99481481e-01, 7.70370370e-03, 3.00000000e-03, 2.89259259e-02, 2.30370370e-02, 5.92592593e-04 },
        { 8.53333333e-02, 6.30666667e-01, 2.82666667e-01, 1.33333333e-03, 0.00000000e+00, 0.00000000e+00, 
          0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00 },
        { 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 5.55555556e-02, 0.00000000e+00, 0.00000000e+00, 
          0.00000000e+00, 7.22222222e-01, 0.00000000e+00, 0.00000000e+00, 0.00000000e+00, 2.22222222e-01 }
    };

    // Expected tensor product values for three evaluation points derived with respect to r: dNi/dr(r) * Nj(s)
    std::vector<std::vector<double>> expectedBasis_10 =
    {
        { -6.00000000e-02, -5.78518519e-01, -4.60740741e-01, -1.18518519e-02,  3.00000000e-02,  2.89259259e-01,
           2.30370370e-01,  5.92592593e-03,  3.00000000e-02,  2.89259259e-01,  2.30370370e-01,  5.92592593e-03 },
        { -2.84444444e-01, -2.10222222e+00, -9.42222222e-01, -4.44444444e-03,  2.84444444e-01,  2.10222222e+00, 
           9.42222222e-01,  4.44444444e-03,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00 },
        { -0.00000000e+00, -0.00000000e+00, -0.00000000e+00, -5.55555556e-01, -0.00000000e+00, -0.00000000e+00, 
          -0.00000000e+00, -5.55555556e-01,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  1.11111111e+00 }
    };

    // Expected tensor product values for three evaluation points derived with respect to s: Ni(r) * dNj/ds(s)
    std::vector<std::vector<double>> expectedBasis_01 =
    {
        { -8.40000000e-02, -1.46222222e-01,  2.05333333e-01,  2.48888889e-02, -2.73000000e-01, -4.75222222e-01,
           6.67333333e-01,  8.08888889e-02, -2.10000000e-02, -3.65555556e-02,  5.13333333e-02,  6.22222222e-03 },
        { -4.48000000e-01, -4.76000000e-01,  8.96000000e-01,  2.80000000e-02, -0.00000000e+00, -0.00000000e+00, 
           0.00000000e+00,  0.00000000e+00, -0.00000000e+00, -0.00000000e+00,  0.00000000e+00,  0.00000000e+00 },
        {  0.00000000e+00,  0.00000000e+00, -2.33333333e-01,  2.33333333e-01,  0.00000000e+00,  0.00000000e+00, 
          -3.03333333e+00,  3.03333333e+00,  0.00000000e+00,  0.00000000e+00, -9.33333333e-01,  9.33333333e-01 }
    };
    
    std::array<size_t, 2> numberOfElements { 5, 7 };
    std::array<size_t, 2> polynomialDegrees { 2, 3 };
    std::array<size_t, 2> continuities { 1, 2 };

    std::array<double, 2> lengths { 3.0, 5.0 };
    std::array<double, 2> origin { -2.0, 4.0 };
    
    auto dummy_integrator = []( size_t ) { return IntegrationPoints { }; };

    auto mesh = BSplineFiniteElementPatch( numberOfElements, polynomialDegrees, continuities, 
                                           lengths, origin, dummy_integrator );

    std::vector<std::vector<double>> computedBasis_00( 3 );
    std::vector<std::vector<double>> computedBasis_10( 3 );
    std::vector<std::vector<double>> computedBasis_01( 3 );

    // Evaluation point 1
    REQUIRE_NOTHROW( computedBasis_00[0] = mesh.evaluateActiveBasisAt( { 0.0,  5.0 }, { 0, 0 } ) );
    REQUIRE_NOTHROW( computedBasis_01[0] = mesh.evaluateActiveBasisAt( { 0.0,  5.0 }, { 0, 1 } ) );
    REQUIRE_NOTHROW( computedBasis_10[0] = mesh.evaluateActiveBasisAt( { 0.0,  5.0 }, { 1, 0 } ) );

    // Evaluation point 2
    REQUIRE_NOTHROW( computedBasis_00[1] = mesh.evaluateActiveBasisAt( { -2.0,  7.0 }, { 0, 0 } ) );
    REQUIRE_NOTHROW( computedBasis_01[1] = mesh.evaluateActiveBasisAt( { -2.0,  7.0 }, { 0, 1 } ) );
    REQUIRE_NOTHROW( computedBasis_10[1] = mesh.evaluateActiveBasisAt( { -2.0,  7.0 }, { 1, 0 } ) );

    // Evaluation point 3
    REQUIRE_NOTHROW( computedBasis_00[2] = mesh.evaluateActiveBasisAt( { -1.0, 9.0 }, { 0, 0 } ) );
    REQUIRE_NOTHROW( computedBasis_01[2] = mesh.evaluateActiveBasisAt( { -1.0, 9.0 }, { 0, 1 } ) );
    REQUIRE_NOTHROW( computedBasis_10[2] = mesh.evaluateActiveBasisAt( { -1.0, 9.0 }, { 1, 0 } ) );
    
    size_t expectedSize = 12; // (p0 + 1) * (p1 + 1)

    for( size_t iPoint = 0; iPoint < 3; ++iPoint )
    {
        REQUIRE( computedBasis_00[iPoint].size( ) == expectedSize );
        REQUIRE( computedBasis_01[iPoint].size( ) == expectedSize );
        REQUIRE( computedBasis_10[iPoint].size( ) == expectedSize );

        for( size_t iFunction = 0; iFunction < expectedSize; ++iFunction )
        {
            CHECK( computedBasis_00[iPoint][iFunction] == Approx( expectedBasis_00[iPoint][iFunction] ).margin( 1e-12 ) );
            CHECK( computedBasis_01[iPoint][iFunction] == Approx( expectedBasis_01[iPoint][iFunction] ).margin( 1e-12 ) );
            CHECK( computedBasis_10[iPoint][iFunction] == Approx( expectedBasis_10[iPoint][iFunction] ).margin( 1e-12 ) );

        } // iFunction
    } // iPoint

} // BSplineFiniteElementPatch_evaluateActiveBasisAt_test



TEST_CASE( "BSplineFiniteElementPatch_integrateElementSystem_test" )
{
    auto gaussIntegrator = []( size_t size ) -> IntegrationPoints
    {
        if( size == 3 )
        {
            return { { { -0.77459667, 0.0,        0.77459667 },
                       {  0.55555556, 0.88888889, 0.55555556 } } };
        }
        if( size == 4 )
        {
            return { { { -0.86113631, -0.33998104,  0.33998104,  0.86113631 },
                       {  0.34785485,  0.65214515,  0.65214515,  0.34785485 } } };
        }
        else
        {
            throw std::runtime_error( "Unexpected integration order." );
        }
    };

    auto testSourceFunction = []( double x, double y )
    {
        return x * y;
    };

    std::vector<std::vector<double>> expectedElementMatrix = 
    {
        {  8.26808398e-03,  1.61542329e-02,  7.86167790e-04, -4.06897207e-04,  8.46595810e-03, -2.47711646e-03, 
          -1.71930840e-02, -1.19655140e-03, -9.84041956e-04, -6.67711643e-03, -4.59308391e-03, -1.46551399e-04 },
        {  1.61542329e-02,  9.37492510e-02,  6.30478840e-02,  6.59744262e-04, -2.47711646e-03, -2.02746255e-02, 
          -5.29906089e-02, -1.10632055e-02, -6.67711643e-03, -4.02246253e-02, -3.68906085e-02, -3.01320547e-03 },
        {  7.86167790e-04,  6.30478840e-02,  1.05437416e-01,  1.26068406e-02, -1.71930840e-02, -5.29906089e-02, 
          -2.09853745e-02,  2.29913046e-04, -4.59308391e-03, -3.68906085e-02, -4.47853743e-02, -4.67008694e-03 },
        { -4.06897207e-04,  6.59744262e-04,  1.26068406e-02,  3.67470399e-03, -1.19655140e-03, -1.10632055e-02, 
           2.29913046e-04,  3.76264804e-03, -1.46551399e-04, -3.01320547e-03, -4.67008694e-03, -4.37351981e-04 },
        {  8.46595810e-03, -2.47711646e-03, -1.71930840e-02, -1.19655140e-03,  4.60680842e-02,  3.29542330e-02, 
          -4.96138324e-02, -4.60689722e-03,  8.46595810e-03, -2.47711646e-03, -1.71930840e-02, -1.19655140e-03 },
        { -2.47711646e-03, -2.02746255e-02, -5.29906089e-02, -1.10632055e-02,  3.29542330e-02,  1.73549252e-01,
          -1.35211593e-03, -3.15402558e-02, -2.47711646e-03, -2.02746255e-02, -5.29906089e-02, -1.10632055e-02 },
        { -1.71930840e-02, -5.29906089e-02, -2.09853745e-02,  2.29913046e-04, -4.96138324e-02, -1.35211593e-03, 
           2.00637416e-01,  3.22068407e-02, -1.71930840e-02, -5.29906089e-02, -2.09853745e-02,  2.29913046e-04 },
        { -1.19655140e-03, -1.10632055e-02,  2.29913046e-04,  3.76264804e-03, -4.60689722e-03, -3.15402558e-02, 
           3.22068407e-02,  2.04747041e-02, -1.19655140e-03, -1.10632055e-02,  2.29913046e-04,  3.76264804e-03 },
        { -9.84041956e-04, -6.67711643e-03, -4.59308391e-03, -1.46551399e-04,  8.46595810e-03, -2.47711646e-03, 
          -1.71930840e-02, -1.19655140e-03,  8.26808398e-03,  1.61542329e-02,  7.86167790e-04, -4.06897207e-04 },
        { -6.67711643e-03, -4.02246253e-02, -3.68906085e-02, -3.01320547e-03, -2.47711646e-03, -2.02746255e-02,
          -5.29906089e-02, -1.10632055e-02,  1.61542329e-02,  9.37492510e-02,  6.30478840e-02,  6.59744262e-04 },
        { -4.59308391e-03, -3.68906085e-02, -4.47853743e-02, -4.67008694e-03, -1.71930840e-02, -5.29906089e-02,
          -2.09853745e-02,  2.29913046e-04,  7.86167790e-04,  6.30478840e-02,  1.05437416e-01,  1.26068406e-02 },
        { -1.46551399e-04, -3.01320547e-03, -4.67008694e-03, -4.37351981e-04, -1.19655140e-03, -1.10632055e-02,
           2.29913046e-04,  3.76264804e-03, -4.06897207e-04,  6.59744262e-04,  1.26068406e-02,  3.67470399e-03 }
    };

    std::vector<double> expectedElementVector = 
    {
        { -2.71045920e-02, -1.95578233e-01, -2.10459185e-01, -1.96641158e-02, -9.54081637e-02, -6.88435377e-01,
          -7.40816330e-01, -6.92176874e-02, -2.05994899e-02, -1.48639457e-01, -1.59948980e-01, -1.49447280e-02 }
    };

    auto mesh = BSplineFiniteElementPatch( { 5, 7 }, { 2, 3 }, { 1, 2 }, { 3.0, 5.0 }, { -2.0, 4.0 }, gaussIntegrator );

    auto computedElementSystem = mesh.integrateElementSystem( { 1, 1 }, testSourceFunction );

    // Check the element matrix
    const auto& computedElementMatrix = std::get<0>( computedElementSystem );

    REQUIRE( computedElementMatrix.size1( ) == expectedElementMatrix.size( ) );
    REQUIRE( computedElementMatrix.size2( ) == expectedElementMatrix.size( ) );
    
    for( size_t i = 0; i < expectedElementMatrix.size( ); ++i )
    {
        for( size_t j = 0; j < expectedElementMatrix.size( ); ++j )
        {
            CHECK( computedElementMatrix( i, j ) == Approx( expectedElementMatrix[i][j] ) );
        }
    }

    // Check the element source vector
    const auto& computedElementRhs = std::get<1>( computedElementSystem );

    REQUIRE( computedElementRhs.size( ) == expectedElementVector.size( ) );

    for( size_t i = 0; i < expectedElementVector.size( ); ++i )
    {
        CHECK( computedElementRhs[i] == Approx( expectedElementVector[i] ) );
    }
}


} // namespace splinekernel
} // namespace cie
