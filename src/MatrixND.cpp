#include "MatrixND.h"

//--Starting point for methods of struct OperatingDimensions_t--
OperatingDimensions_t::OperatingDimensions_t(void)
{
	da = 1;
	db = 2;
}

OperatingDimensions_t::OperatingDimensions_t(UINT16 da, UINT16 db)
{
	set(da, db);
}

void OperatingDimensions_t::set(UINT16 da, UINT16 db)
{
	{
		//equivalent values will not change the values as this is not a valid case
		if (da == db) return;
		if (da < db)
		{
			this->da = da;
			this->db = db;
			return;
		}
		this->da = db;
		this->db = da;
		return;
	}
}

//---------Starting point for methods of class MatrixND----------
MatrixND::MatrixND(std::vector<UINT32> dimensions)
{
	m_iDimensionality = dimensions.size();
	m_iElements = 1;
	m_piDimensions = new UINT32[m_iDimensionality];
	for (UINT16 i = 0; i < m_iDimensionality; i++)
	{
		m_piDimensions[i] = dimensions.at(i);
		m_iElements *= m_piDimensions[i];
	}
	m_pfData = new float[m_iElements];
	for (UINT32 j = 0; j < m_iElements; j++)
	{
		m_pfData[j] = 0.0f;
	}
	m_modPrevent = 0;
}

MatrixND::~MatrixND(void)
{
}

float& MatrixND::at(UINT32 index)
{
	if (isInMatrix(index))
	{
		return m_pfData[index];
	}
	return m_modPrevent;
}

float& MatrixND::at(const std::vector<UINT32>& position)
{
	if (isInMatrix(position))
	{
		return m_pfData[getIndexFromPosition(position)];
	}
	return m_modPrevent;
}

float& MatrixND::atFast(UINT32 index)
{
	return m_pfData[index];
}

float& MatrixND::atFast(UINT32* position)
{
	return m_pfData[getIndexFromPositionFast(position)];
}

MatrixND MatrixND::transpose(MatrixND matIn, OperatingDimensions_t dims)
{
	if (!matIn.isInMatrix(dims.da) || !matIn.isInMatrix(dims.db)) 
		return matIn;
	std::vector<UINT32> dimensions(matIn.m_iDimensionality);
	for (UINT32 i = 0; i < matIn.m_iDimensionality; i++)
	{
		dimensions.at(i) = matIn.m_piDimensions[i];
	}
	dimensions.at(dims.da - 1) = matIn.m_piDimensions[dims.db - 1];
	dimensions.at(dims.db - 1) = matIn.m_piDimensions[dims.da - 1];
	MatrixND matOut(dimensions);
	for (UINT32 j = 0; j < matOut.m_iElements; j++)
	{
		UINT32* position = matOut.getPositionFromIndexFast(j);
		UINT32* newPosition = matOut.getPositionFromIndexFast(j);
		//Swaps the position values at the appropriate dimensions
		newPosition[dims.da - 1] = position[dims.db - 1];
		newPosition[dims.db - 1] = position[dims.da - 1];
		matOut.atFast(j) = matIn.atFast(newPosition);
	}
	return matOut;
}

MatrixND MatrixND::generateIdentity(std::vector<UINT32>& dimensions, OperatingDimensions_t dims)
{
	MatrixND identity(dimensions);
	if (dimensions.at(dims.da - 1) == dimensions.at(dims.db - 1))
	{
		for (UINT32 i = 0; i < identity.m_iElements; i++)
		{
			UINT32* position = identity.getPositionFromIndexFast(i);
			if (position[dims.da - 1] == position[dims.db - 1])
			{
				identity.atFast(i) = 1.0f;
			}
		}
	}
	return identity;
}

//---------------------------Operations--------------------------

MatrixND& MatrixND::scalarMultiply(float multiple)
{
	for (UINT32 i = 0; i < m_iElements; i++)
	{
		m_pfData[i] *= multiple;
	}
	return *this;
}

MatrixND& MatrixND::add(MatrixND other)
{
	if (compareDimensions(other))
	{
		for (UINT32 i = 0; i < m_iElements; i++)
		{
			m_pfData[i] += other.m_pfData[i];
		}
	}
	return *this;
}

MatrixND& MatrixND::subtract(MatrixND other)
{
	if (compareDimensions(other))
	{
		for (UINT32 i = 0; i < m_iElements; i++)
		{
			m_pfData[i] -= other.m_pfData[i];
		}
	}
	return *this;
}

MatrixND& MatrixND::multiply(MatrixND other)
{
	if (multipliable(other))
	{
		UINT32 n = m_piDimensions[m_OperatingDimensions.da];
		std::vector<UINT32> dimensions(this->m_iDimensionality);
		for (UINT32 dimension = 0; dimension < m_iDimensionality; dimension++)
		{
			dimensions.at(dimension) = this->m_piDimensions[dimension];
		}
		dimensions.at(m_OperatingDimensions.da - 1) = this->m_piDimensions[m_OperatingDimensions.da - 1];
		dimensions.at(m_OperatingDimensions.db - 1) = this->m_piDimensions[m_OperatingDimensions.db - 1];

		MatrixND matOut(dimensions);
		for (UINT32 i = 0; i < matOut.m_iElements; i++)
		{
			float sum = 0;
			for (UINT32 x = 0; x < n; x++)
			{
				UINT32* positionA = getPositionFromIndexFast(i);
				UINT32* positionB = getPositionFromIndexFast(i);
				positionA[m_OperatingDimensions.db - 1] = x + 1;
				positionB[m_OperatingDimensions.da - 1] = x + 1;
				sum += (this->atFast(positionA)) * (other.atFast(positionB));
			}
			matOut.atFast(i) = sum;
		}
		matOut.copy(this);
	}
	return *this;
}

MatrixND& MatrixND::outerProduct(MatrixND other)
{
	UINT32 dimensionality = this->m_iDimensionality + other.m_iDimensionality;
	std::vector<UINT32> dimensions(dimensionality);
	for (UINT32 i = 0; i < dimensionality; i++)
	{
		if (i < this->m_iDimensionality)
		{
			dimensions.at(i) = this->m_piDimensions[i];
		}
		else
		{
			dimensions.at(i) = other.m_piDimensions[i - this->m_iDimensionality];
		}
	}

	MatrixND matOut(dimensions);
	for (UINT32 j = 0; j < this->m_iElements; j++)
	{
		UINT32* firstPosition = this->getPositionFromIndexFast(j);
		float firstValue = this->atFast(j);

		for (UINT32 k = 0; k < other.m_iElements; k++)
		{
			UINT32* secondPosition = this->getPositionFromIndexFast(k);
			float secondValue = other.atFast(k);

			UINT32* resultingPosition = new UINT32[dimensionality];
			//Appending the two arrays to get the resulting position
			for (UINT32 index = 0; index < dimensionality; index++)
			{
				if (index < this->m_iDimensionality)
				{
					resultingPosition[index] = firstPosition[index];
				}
				else
				{
					resultingPosition[index] = secondPosition[index - this->m_iDimensionality];
				}

			}
			matOut.atFast(resultingPosition) = firstValue * secondValue;
		}
	}
	matOut.copy(this);
	return *this;
}

bool MatrixND::equals(MatrixND other) const
{
	if (!compareDimensions(other))
		return false;
	if (this->m_iElements != other.m_iElements)
		return false;
	for (UINT32 i = 0; i < m_iElements; i++)
	{
		if (this->m_pfData[i] != other.m_pfData[i])
			return false;
	}
	return true;
}

//-------------------------Operators------------------------------

MatrixND operator^(MatrixND mat1, OperatingDimensions_t dims)
{
	return MatrixND::transpose(mat1, dims);
}

MatrixND& MatrixND::operator+=(MatrixND other)
{
	return add(other);
}

MatrixND& MatrixND::operator-=(MatrixND other)
{
	return subtract(other);
}

MatrixND& MatrixND::operator*=(float multiple)
{
	return scalarMultiply(multiple);
}

MatrixND MatrixND::operator*=(MatrixND other)
{
	return multiply(other);
}

MatrixND& operator+(MatrixND mat1, MatrixND mat2)
{
	return mat1.add(mat2);
}

MatrixND& operator-(MatrixND mat1, MatrixND mat2)
{
	return mat1.subtract(mat2);
}

MatrixND& operator*(float& multiple, MatrixND mat)
{
	return mat.scalarMultiply(multiple);
}

MatrixND& operator*(MatrixND mat, float& multiple)
{
	return mat.scalarMultiply(multiple);
}

MatrixND operator*(MatrixND mat1,MatrixND mat2)
{
	return mat1.multiply(mat2);
}

bool operator==(MatrixND mat1, MatrixND mat2)
{
	return mat1.equals(mat2);
}

bool operator!=(MatrixND mat1, MatrixND mat2)
{
	return !mat1.equals(mat2);
}

//------------------------Utilities----------------------------

void MatrixND::copy(MatrixND* target)
{
	target->m_iDimensionality = this->m_iDimensionality;
	target->m_iElements = this->m_iElements;
	target->m_piDimensions = NULL;
	target->m_pfData = NULL;
	target->m_piDimensions = new UINT32[target->m_iDimensionality];
	target->m_pfData = new float[target->m_iElements];
	memcpy(target->m_piDimensions, this->m_piDimensions, sizeof(UINT32) * this->m_iDimensionality);
	memcpy(target->m_pfData, this->m_pfData, sizeof(UINT32) * this->m_iElements);
}

void MatrixND::setOperatingDimensions(UINT16 da, UINT16 db)
{
	m_OperatingDimensions.set(da, db);
}

//-------------------------Positioners-----------------------------

std::vector<UINT32> MatrixND::getPositionFromIndex(UINT32 index) const
{
	std::vector<UINT32> position(m_iDimensionality);
	UINT32 product = 1;
	UINT32 tempIndex = index;
	for (UINT16 i = 0; i < m_iDimensionality; i++)
	{
		product *= m_piDimensions[i];
		position.at(i) = 1;
	}
	UINT32 currentDimension = m_iDimensionality - 1;
	while (currentDimension < m_iDimensionality)
	{
		product /= m_piDimensions[currentDimension];
		//We check if the value is less than the current
		//In essence checking for integer wrapping when 
		//the value tries to be negative in unsigned type
		while ((tempIndex - product) < tempIndex)
		{
			position.at(currentDimension) += 1;
			tempIndex -= product;
		}
		if (tempIndex == 0) break;
		currentDimension--;
	}
	return position;
}

UINT32 MatrixND::getIndexFromPosition(std::vector<UINT32> pos) const
{
		UINT32 index = 0;
		UINT32 product = 1;
		for (UINT16 i = 0; i < m_iDimensionality; i++)
		{
			product *= m_piDimensions[i];
		}
		for (UINT16 j = m_iDimensionality - 1; j <= m_iDimensionality - 1; j--)
		{
			product /= m_piDimensions[j];
			index += (pos.at(j) - 1) * product;
		}
		return index;
}

UINT32* MatrixND::getPositionFromIndexFast(UINT32 index) const
{
	UINT32* position = new UINT32[m_iDimensionality];
	UINT32 product = 1;
	UINT32 tempIndex = index;
	for (UINT16 i = 0; i < m_iDimensionality; i++)
	{
		product *= m_piDimensions[i];
		position[i] = 1;
	}
	UINT32 currentDimension = m_iDimensionality - 1;
	while (currentDimension < m_iDimensionality)
	{
		product /= m_piDimensions[currentDimension];
		//We check if the value is less than the current
		//In essence checking for integer wrapping when 
		//the value tries to be negative in unsigned type
		while ((tempIndex - product) < tempIndex)
		{
			position[currentDimension] += 1;
			tempIndex -= product;
		}
		if (tempIndex == 0) break;
		currentDimension--;
	}
	return position;
}

UINT32 MatrixND::getIndexFromPositionFast(UINT32* pos) const
{
	UINT32 index = 0;
	UINT32 product = 1;
	for (UINT16 i = 0; i < m_iDimensionality; i++)
	{
		product *= m_piDimensions[i];
	}
	for (UINT16 j = m_iDimensionality - 1; j <= m_iDimensionality - 1; j--)
	{
		product /= m_piDimensions[j];
		index += (pos[j] - 1) * product;
	}
	return index;
}

//------------------------General Checkers------------------------

bool MatrixND::isInMatrix(std::vector<UINT32> pos) const
{
	if (pos.size() != this->m_iDimensionality)
		return false;
	for (UINT16 i = 0; i < m_iDimensionality; i++)
	{
		if (pos.at(i) > m_piDimensions[i])
			return false;
	}
	return true;
}

bool MatrixND::isInMatrix(UINT32 index) const
{
		return index <= m_iElements;
}

bool MatrixND::dimensionExists(const UINT16& dimension) const
{
	return dimension <= m_iDimensionality;
}

bool MatrixND::compareDimensions(MatrixND other) const
{
	if (m_iDimensionality != other.m_iDimensionality)
		return false;
	for (UINT16 i = 0; i < m_iDimensionality; i++)
	{
		if (m_piDimensions[i] != other.m_piDimensions[i])
			return false;
	}
	return true;
}

bool MatrixND::multipliable(MatrixND other) const
{
	if (m_iDimensionality != other.m_iDimensionality)
		return false;
	if (this->m_piDimensions[m_OperatingDimensions.db - 1] != other.m_piDimensions[m_OperatingDimensions.da - 1])
		return false;
	for (UINT16 d = 0; d < m_iDimensionality; d++)
	{
		if (d == m_OperatingDimensions.da - 1 || d == m_OperatingDimensions.db - 1)
			continue;
		if (this->m_piDimensions[d] != other.m_piDimensions[d])
			return false;
	}
	return true;
}