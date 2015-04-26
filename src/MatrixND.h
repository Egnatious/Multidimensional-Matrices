/*****************************************Comment**********************************************
*Header file for MatrixND
*Purpose:  To be able to represent and operate on Multidimensional Matrices
*Author(s): Egnatious (Jordan Ericksen)
*Date Created: 3/30/15 
*Date Last Modified: 4/23/15
*Version: 1.0
****************************************End Comment********************************************/
#pragma once

#define DllExport   __declspec( dllexport )

#include<vector>
//defined to prevent dependency on intsafe.h for Mac and Linux platforms
typedef unsigned int UINT32;
typedef unsigned short UINT16;

/*This struct defines the dimensions the matrix will do its calculations in.
By default these values are one and two. Although the set function will
automatically place the values in proper order the lowest value should be first.*/
struct OperatingDimensions_t
{
	UINT16 da;
	UINT16 db;

	DllExport OperatingDimensions_t(void);
	DllExport OperatingDimensions_t(UINT16 da, UINT16 db);

	DllExport void set(UINT16 da, UINT16 db);
};

class MatrixND
{
public:
	//Constructors
	DllExport MatrixND(std::vector<UINT32> dimensions);
	DllExport ~MatrixND(void);
private:
	//Class Members
	float* m_pfData;
	UINT32* m_piDimensions;
	UINT16 m_iDimensionality;
	UINT32 m_iElements;
	OperatingDimensions_t m_OperatingDimensions;
	/*This is to keep someone from modifying a non-existent
	reference and to maintain external memory security*/
	float m_modPrevent;
public:
	//Public functions
	DllExport float& at(UINT32 index);
	DllExport float& at(const std::vector<UINT32>& position);

	/*Does not check whether the values are in matrix and trusts the programmer
	Useful for performance in trusted code i.e. The multiplication function*/
	DllExport inline float& atFast(UINT32 index);
	DllExport inline float& atFast(UINT32* position);

	DllExport static MatrixND generateIdentity(std::vector<UINT32>& dimensions, OperatingDimensions_t dims);
	DllExport static MatrixND transpose(MatrixND matIn, OperatingDimensions_t dims);

	DllExport inline MatrixND& scalarMultiply(float multiple);
	DllExport inline MatrixND& add(MatrixND other);
	DllExport inline MatrixND& subtract(MatrixND other);

	DllExport MatrixND& multiply(MatrixND other);

	DllExport bool equals(MatrixND other) const;

	DllExport MatrixND& outerProduct(MatrixND other);

	DllExport MatrixND& operator+=(MatrixND other);
	DllExport MatrixND& operator-=(MatrixND other);
	DllExport MatrixND& operator*=(float multiple);
	DllExport MatrixND operator*=(MatrixND other);

	//Using the bitwise xor to signify transpose with operators
	DllExport friend MatrixND operator^(MatrixND mat1, OperatingDimensions_t dims);

	DllExport friend MatrixND& operator+(MatrixND mat1, MatrixND mat2);
	DllExport friend MatrixND& operator-(MatrixND mat1, MatrixND mat2);
	DllExport friend MatrixND& operator*(float& multiple, MatrixND mat);
	DllExport friend MatrixND& operator*(MatrixND mat, float* multiple);
	DllExport friend MatrixND operator*(MatrixND mat1, MatrixND mat2);

	DllExport friend bool operator==(MatrixND mat1, MatrixND mat2);
	DllExport friend bool operator!=(MatrixND mat1, MatrixND mat2);

	DllExport void copy(MatrixND* target);
	DllExport void setOperatingDimensions(UINT16 da, UINT16 db);

	//Functions only appears in header
	DllExport inline UINT16 getDimensionality(void) const{return m_iDimensionality;}
	DllExport inline UINT32 getElements(void) const{return m_iElements;}
	DllExport inline UINT32* const getDimensions(void) const{return m_piDimensions;}
private:
	//Private Functions
	std::vector<UINT32> getPositionFromIndex(UINT32 index) const;
	UINT32 getIndexFromPosition(std::vector<UINT32> pos) const;
	UINT32* getPositionFromIndexFast(UINT32 index) const;
	UINT32 getIndexFromPositionFast(UINT32* pos) const;

	inline bool isInMatrix(std::vector<UINT32> pos) const;
	inline bool isInMatrix(UINT32 index) const;
	inline bool dimensionExists(const UINT16& dimension) const;
	inline bool compareDimensions(MatrixND other) const;

	bool multipliable(MatrixND other) const;
};

/***********************************************Comment*********************************************************
*The following links will show the papers used to define the rules being used in the program
*
*Multidimensional Matrix Mathematics:
*	Part 1: Notation, Representation, and Simplification
*		http://www.iaeng.org/publication/WCE2010/WCE2010_pp1824-1828.pdf
*	Part 2: Multidimensional Matrix Equality, Addition, Subtraction, and Multiplication
*		http://www.iaeng.org/publication/WCE2010/WCE2010_pp1829-1833.pdf
*	Part 3: Multidimensional Null and Identity Matrices and Multidimensional Matrix Outer and Inner Products
*		http://www.iaeng.org/publication/WCE2010/WCE2010_pp1834-1837.pdf
*	Part 4: Multidimensional Matrix Transpose, Symmetry, Antisymmetry, Determinant, and Inverse
*		http://www.iaeng.org/publication/WCE2010/WCE2010_pp1838-1841.pdf
*	Part 5: Algebraic Laws
*		http://www.iaeng.org/publication/WCE2010/WCE2010_pp1842-1847.pdf
*	Part 6: Solving Systems of Linear Equations and Multidimensional Matrix Calculus
*		http://www.iaeng.org/publication/WCE2010/WCE2010_pp1848-1850.pdf
***********************************************End Comment*******************************************************/