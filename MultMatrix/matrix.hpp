#ifndef LR3_MATRIX_H
#define LR3_MATRIX_H

#include <vector>
#include "types.hpp"

struct Matrix {
private:
    u32 m, n;

public:
	enum EMatrixErrors {
		ME_None			= 0x0000,
		ME_Init			= 0x0001,
		ME_Dimensions	= 0x0002,
		ME_NotSquare	= 0x0004,
		ME_Singular		= 0x0008,

		ME_Size			= ME_Dimensions
	};

    u32 error = ME_None;

	void ClearError( void ) {
		error = ME_None;
	}

	// Ultrafast, but you must be a bit more careful to use it. C-style, in other words.
	double **data = NULL;

	// All the pointers are saved as-is (without copying the content).
    Matrix( u32 x, u32 y );
    Matrix( u32 x, u32 y, double **cells );
    Matrix( u32 x, u32 y, float **cells );
	Matrix( const Matrix& ); // Copying.

     explicit Matrix( u32 squareside ) : Matrix( squareside, squareside ) {};
    Matrix( u32 squareside, double **cells ) : Matrix( squareside, squareside, cells ) {};

	//Matrix( std::initializer_list<std::initializer_list<double>> );
	explicit Matrix( std::vector<std::vector<double>> );

	Matrix& operator= ( const Matrix& ) noexcept;// Moving.
	//Matrix* operator=(Matrix&) noexcept;
	Matrix operator= ( const std::vector<std::vector<double>>& );
    Matrix operator+ ( const Matrix& );
    Matrix operator- ( const Matrix& );
    Matrix operator* ( const Matrix& );
    Matrix operator* ( double );

	const Matrix* operator* () { return this; }

	bool operator== ( const Matrix& );
	bool operator!= ( const Matrix& );

    u32 GetRows( void ) { return n; }
    u32 GetCols( void ) { return m; }

	void DebugPrint( void );

/*
	double Det( void );
	void Transpose( void );
	void Inverse( void );

	void Rotate( double angle, double pitch );
	void Move( double x, double y, double z );
	//void Project( ? );
*/

	~Matrix();

    static Matrix MoveBy(float dx, float dy, float dz);

    // angle in rad
    static Matrix RotateX(float angle);
    static Matrix RotateY(float angle);
    static Matrix RotateZ(float angle);

};

#endif // of #ifndef LR3_MATRIX_H
