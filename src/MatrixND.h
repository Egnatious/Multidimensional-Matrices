/*****************************************Comment**********************************************
*Header file for MatrixND
*Purpose:  To be able to represent and operate on Multidimensional Matrices
*Author(s): Egnatious (Jordan Ericksen)
*Date Created: 3/30/15 
*Date Last Modified: 4/6/15
*Version: 0.1
*Notes: This class uses column ordered matrices because that was what was used in the papers
*linked on the bottom of the header file
****************************************End Comment********************************************/
#pragma once

#ifdef _WIN32
#include <intsafe.h>
#elif defined(__linux__) || (defined(__APPLE__) && defined(__MACH__))
typedef UINT16 unsigned char;
typedef UINT32 unsigned int;
#define UINT32_MAX      0xffffffffui32
#endif
#include <vector>

using namespace std;


class MatrixND
{
public:
	MatrixND(vector<UINT32> dimensions);
	~MatrixND(void);
	
	/*This struct defines the dimensions the matrix will do its calculations in.
	By default these values are one and two. Although the set function will
	automatically place the values in proper order the lowest value should be first.*/
	struct OperatingDimensions_t
	{
		UINT16 da;
		UINT16 db;

		OperatingDimensions_t(void);
		OperatingDimensions_t(UINT16 da, UINT16 db);

		void set(UINT16 da, UINT16 db);
	};
private:
	float* m_pfData;
	UINT32* m_piDimensions;
	UINT16 m_iDimensionality;
	UINT32 m_iElements;
	OperatingDimensions_t m_OperatingDimensions;
	/*This is to keep someone from modifying a non-existent
	reference and maintain external memory security*/
	float m_modPrevent;
public:
	inline float& at(UINT32 index);
	inline float& at(vector<UINT32> position);

	inline void setOperatingDimensions(UINT16 da, UINT16 db);

	static MatrixND generateIdentity(vector<UINT32>& dimensions, OperatingDimensions_t dims);
	static MatrixND transpose(MatrixND matIn, OperatingDimensions_t dims);

	inline MatrixND& scalarMultiply(float multiple);
	inline MatrixND& add(MatrixND other);
	inline MatrixND& subtract(MatrixND other);

	MatrixND multiply(MatrixND other);

	MatrixND outerProduct(MatrixND other);
	MatrixND& operator+=(MatrixND other);
	MatrixND& operator-=(MatrixND other);
	MatrixND& operator*=(float multiple);
	MatrixND operator*=(MatrixND other);

	friend MatrixND& operator+(MatrixND mat1, MatrixND mat2);
	friend MatrixND& operator-(MatrixND mat1, MatrixND mat2);
	friend MatrixND& operator*(float& multiple, MatrixND mat);
	friend MatrixND& operator*(MatrixND mat, float* multiple);
	friend MatrixND operator*(MatrixND mat1, MatrixND mat2);
	//Using the bitwise xor to signify transpose with operators
	friend MatrixND operator^(MatrixND mat1, OperatingDimensions_t dims);

	inline UINT32 getElements(void){return m_iElements;}
private:
	vector<UINT32> getPositionFromIndex(UINT32 index);
	UINT32 getIndexFromPosition(vector<UINT32> pos);
	
	bool multipliable(MatrixND other) const;

	inline bool isInMatrix(vector<UINT32> pos) const;
	inline bool isInMatrix(UINT32 index) const;
	inline bool dimensionExists(const UINT16& dimension) const;
	inline bool compareDimensions(MatrixND other) const;
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