#include <cstdio>
#include "matrix.hpp"
#include "types.hpp"
#include <cmath>
#include <iostream>

Matrix::~Matrix() {
    for ( u32 i = 0; i < n; i++ )
		delete[] data[ i ];

    delete[] data;
}


Matrix::Matrix( const Matrix& src ) {
	n = src.n;
	m = src.m;

	data = new double*[ n ];

    for ( u32 i = 0; i < n; i++ )
		data[ i ] = new double[ m ];

    for ( u32 i = 0; i < n; i++ )
        for ( u32 j = 0; j < m; j++ )
			data[ i ][ j ] = src.data[ i ][ j ];
}


Matrix::Matrix( u32 x, u32 y ) {
	data = new double*[ y ];

    for ( u32 i = 0; i < y; i++ )
		data[ i ] = new double[ x ];

	m = x;
	n = y;

	//printf( "Created xy %ux%u: %p. data: %p\n", (unsigned int) x, (unsigned int) y, (void *) this, (void *)data );
}

Matrix::Matrix( u32 x, u32 y, double **cells ) {
	data = cells;
	m = x;
	n = y;
}

Matrix::Matrix( u32 x, u32 y, float **cells ) {
    m = x;
    n = y;
    data = new double*[ n ];

    for ( u32 i = 0; i < n; i++ )
        data[ i ] = new double[ m ];

    for ( u32 i = 0; i < n; i++ )
        for ( u32 j = 0; j < m; j++ )
            data[ i ][ j ] = cells[ i ][ j ];
}

Matrix::Matrix( std::vector<std::vector<double>> cells ) {
	n = cells.size();

	if ( n > 0 ) {
		m = cells[ 0 ].size();

		data = new double*[ n ];

        for ( u32 i = 0; i < n; i++ ) {
			data[ i ] = new double[ m ];

            u32 currowsize = cells[ i ].size();

			if ( currowsize >= m ) {
				// Fill up to m cells.
                for ( u32 j = 0; j < m; j++ )
					data[ i ][ j ] = cells[ i ][ j ];
			} else {
				// Fill all cells and add zeros to the end.
                for ( u32 j = 0; j < currowsize; j++ )
					data[ i ][ j ] = cells[ i ][ j ];

                for ( u32 j = currowsize; j < m; j++ )
					data[ i ][ j ] = 0.0;
			}

		}
	} else {
		m = 0;
		error = ME_Init | ME_Dimensions;
	}

	//printf( "Created A<A<>> %ux%u: %p\n", (unsigned int) m, (unsigned int) n, (void *) this );
} // of Matrix::Matrix( std::vector<std::vector<double>> cells ) {}


// Moving? [McM]: Don't sure how to say to THIS FUCKING STUPID C++ to handle it properly.
Matrix& Matrix::operator= ( const Matrix& right ) noexcept {
	if ( this == &right )
		return *this;

    for ( u32 i = 0; i < n; i++ )
		delete[] data[ i ];

	delete[] data;

	n = right.n;
	m = right.m;

	data = new double*[ n ];

    for ( u32 i = 0; i < n; i++ )
		data[ i ] = new double[ m ];

	for ( u32 i = 0; i < n; i++ ) {
		for ( u32 j = 0; j < m; j++ )
			data[ i ][ j ] = right.data[ i ][ j ];
	}

	//printf( "Assigned %p %ux%u. data: %p\n", (void *) this, (unsigned int) m, (unsigned int) n, (void *)data );
	//DebugPrint();

	return *this;
}

Matrix Matrix::operator= ( const std::vector<std::vector<double>>& right ) {
	return Matrix( right );
}

Matrix Matrix::operator+ ( const Matrix& right ) {
	if ( m != right.m || n != right.n ) {
		error |= ME_Dimensions;
		puts( "Error: dimensions mismatch in sum (A + B)." );
	} else {
		Matrix C( m, n );

        for ( u32 i = 0; i < n; i++ )
            for ( u32 j = 0; j < m; j++ )
				C.data[ i ][ j ] = data[ i ][ j ] + right.data[ i ][ j ];

		return C;
	}

    return Matrix(*this);
}

Matrix Matrix::operator- ( const Matrix& right ) {
    Matrix C( m, n );
	if ( m != right.m || n != right.n ) {
        C.error |= ME_Dimensions;
	} else {

        for ( u32 i = 0; i < n; i++ )
            for ( u32 j = 0; j < m; j++ )
				C.data[ i ][ j ] = data[ i ][ j ] - right.data[ i ][ j ];

	}
    return C;
}

Matrix Matrix::operator* ( const Matrix& right ) {
    Matrix C( right.m, n );

	if ( m != right.n ) {
        C.error |= ME_Dimensions;
		printf( "Error: dimensions mismatch (%lux%lu * %lux%lu).\n", n, m, right.n, right.m );
	} else {

        for ( u32 i = 0; i < n; i++ ) {
            for ( u32 j = 0; j < right.m; j++ ) {
				double sum = 0.0;

                for ( u32 k = 0; k < right.n; k++ )
					sum += data[ i ][ k ] * right.data[ k ][ j ];

				C.data[ i ][ j ] = sum;
			}
		}

	}

    return C;
}

Matrix Matrix::operator* ( double right ) {
    for ( u32 i = 0; i < n; i++ )
        for ( u32 j = 0; j < m; j++ )
			data[ i ][ j ] *= right;

    return Matrix(*this);
}



bool Matrix::operator== ( const Matrix& right ) {
	if ( n != right.n || m != right.m )
		return false;

	for ( u32 i = 0; i < n; i++ )
		for ( u32 j = 0; j < m; j++ )
			if ( data[ i ][ j ] != right.data[ i ][ j ] )
				return false;

	return true;
}

bool Matrix::operator!= ( const Matrix& right ) {
    return !( *this == right );
}

void Matrix::DebugPrint( void ) {
    for ( u32 i = 0; i < n; i++ ) {
		printf( " " );

        for ( u32 j = 0; j < m; j++ )
			std::cout<<data[i][j]<<" ";

		puts( "" ); // A newline.
	}
}

Matrix Matrix::MoveBy(float dx, float dy, float dz)
{
    Matrix res ({
        { 1, 0, 0, 0},
        { 0, 1, 0, 0},
        { 0, 0, 1, 0},
        {dx,dy,dz, 1}
    });
    return res;
}

Matrix Matrix::RotateX(float angle) {
    float sA = sin(angle);
    float cA = cos(angle);

    Matrix res ({
        {1,  0,  0,  0},
        {0, cA, sA,  0},
        {0,-sA, cA,  0},
        {0,  0,  0,  1}
    });
    return res;
}

Matrix Matrix::RotateY(float angle) {
    float sA = sin(angle);
    float cA = cos(angle);

    Matrix res ({
        {cA,  0,-sA,  0},
        { 0,  1,  0,  0},
        {sA,  0, cA,  0},
        { 0,  0,  0,  1}
    });
    return res;
}

Matrix Matrix::RotateZ(float angle) {
    float sA = sin(angle);
    float cA = cos(angle);

    Matrix res ({
        { cA, sA,  0,  0},
        {-sA, cA,  0,  0},
        {  0,  0,  1,  0},
        {  0,  0,  0,  1}
    });
    return res;
}

//Matrix *Matrix::operator=(Matrix& matrix) noexcept {
//
//    if ( this == &matrix )
//        return this;
//
//    std::move(this, matrix)
//
//    for ( u32 i = 0; i < n; i++ )
//        delete[] data[ i ];
//
//    delete[] data;
//    return nullptr;
//}
