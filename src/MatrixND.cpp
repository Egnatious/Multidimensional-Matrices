#include "MatrixND.h"


MatrixND::MatrixND(vector<UINT32> dimensions)
{
	m_iDimensionality = dimensions.size();
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
}

MatrixND::~MatrixND()
{
}

float& MatrixND::at(UINT32 index)
{
	if (isInMatrix(index))
		return m_pfData[index];
	return m_modPrevent;
}

float& MatrixND::at(vector<UINT32> position)
{
	if (isInMatrix(position))
		return m_pfData[getIndexFromPosition(position)];
	return m_modPrevent;
}

void MatrixND::setOperatingDimensions(UINT16 da, UINT16 db)
{
	m_OperatingDimensions.set(da, db);
}

MatrixND MatrixND::transpose(MatrixND matIn, OperatingDimensions_t dims)
{
	if (!matIn.isInMatrix(dims.da) || !matIn.isInMatrix(dims.db)) 
		return matIn;
	vector<UINT32> dimensions(matIn.m_iDimensionality);
	for (UINT32 i = 0; i < matIn.m_iDimensionality; i++)
	{
		dimensions.at(i) = matIn.m_piDimensions[i];
	}
	dimensions.at(dims.da) = matIn.m_piDimensions[dims.db];
	dimensions.at(dims.db) = matIn.m_piDimensions[dims.da];
	MatrixND matOut(dimensions);
	for (UINT j = 0; j < matOut.m_iElements; j++)
	{
		vector<UINT32> position = matOut.getPositionFromIndex(j);
		vector<UINT32> newPosition = matOut.getPositionFromIndex(j);
		//Swaps the position values at the appropriate dimensions
		newPosition.at(dims.da) = position.at(dims.db);
		newPosition.at(dims.db) = position.at(dims.da);
		matOut.at(j) = matIn.at(newPosition);
	}
	return matOut;
}

MatrixND MatrixND::generateIdentity(vector<UINT32>& dimensions, OperatingDimensions_t dims)
{
	MatrixND identity(dimensions);
	if (dimensions.at(dims.da) == dimensions.at(dims.db))
	{
		vector<UINT32> position(identity.m_iDimensionality);
		for (UINT16 k = 0; k < identity.m_iDimensionality; k++) 
			position.at(k) = 0;
		for (UINT16 i = 0; i < dimensions.at(dims.da); i++)
		{
			position.at(dims.da) = i;
			position.at(dims.db) = i;
			identity.at(position) = 1;
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

MatrixND MatrixND::multiply(MatrixND other)
{
	if (multipliable(other))
	{
		UINT32 n = m_piDimensions[m_OperatingDimensions.da];
		vector<UINT32> dimensions(this->m_iDimensionality);
		for (UINT32 dimension = 0; dimension < m_iDimensionality; dimension++)
		{
			dimensions.at(dimension) = this->m_piDimensions[dimension];
		}
		dimensions.at(m_OperatingDimensions.da) = this->m_piDimensions[m_OperatingDimensions.da];
		dimensions.at(m_OperatingDimensions.db) = this->m_piDimensions[m_OperatingDimensions.db];

		MatrixND matOut(dimensions);
		for (UINT32 i = 0; i < matOut.m_iElements; i++)
		{
			float sum = 0;
			for (UINT32 x = 0; x < n; x++)
			{
				vector<UINT32> positionA = getPositionFromIndex(i);
				vector<UINT32> positionB = getPositionFromIndex(i);
				positionA.at(m_OperatingDimensions.db);
				positionB.at(m_OperatingDimensions.da);
				sum += this->at(positionA) * other.at(positionB);
			}
			matOut.at(i) = sum;
		}
		return matOut;
	}
	else
	{
	return *this;
	}
}

MatrixND MatrixND::outerProduct(MatrixND other)
{
	UINT32 dimensionality = this->m_iDimensionality + other.m_iDimensionality;
	vector<UINT32> dimensions(dimensionality);
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
		vector<UINT32> firstPosition = this->getPositionFromIndex(j);
		float firstValue = this->at(j);

		for (UINT32 k = 0; k < other.m_iElements; k++)
		{
			vector<UINT32> secondPosition = other.getPositionFromIndex(k);
			float secondValue = other.at(k);

			vector<UINT32> resultingPosition;
			//Appending the two vectors to get the resulting position
			resultingPosition.insert(resultingPosition.end(), firstPosition.begin(), firstPosition.end());
			resultingPosition.insert(resultingPosition.end(), secondPosition.begin(), secondPosition.end());
			matOut.at(resultingPosition) = firstValue * secondValue;
		}
	}
	return matOut;
}

//------------------------General Checkers------------------------

bool MatrixND::isInMatrix(vector<UINT32> pos) const
{
	if (pos.size() != this->m_iDimensionality)
		return false;
	for (UINT16 i = 0; i < m_iDimensionality; i++)
	{
		if (pos.at(i) != m_piDimensions[i])
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
	if (this->m_piDimensions[m_OperatingDimensions.db] != other.m_piDimensions[m_OperatingDimensions.da])
		return false;
	for (UINT16 d = 0; d < m_iDimensionality; d++)
	{
		if (d == m_OperatingDimensions.da || d == m_OperatingDimensions.db)
			continue;
		if (this->m_piDimensions[d] != other.m_piDimensions[d])
			return false;
	}
	return true;
}

//-------------------------Positioners-----------------------------

vector<UINT32> MatrixND::getPositionFromIndex(UINT32 index)
{
	vector<UINT32> position(m_iDimensionality);
	UINT32 product = 1;
	UINT32 tempIndex = index;
	for (UINT16 i = 0; i < m_iDimensionality - 1; i++)
	{
		product *= m_piDimensions[i];
		position.at(i) = 1;
	}
	UINT32 currentDimension = m_iDimensionality - 1;
	while (currentDimension > 1)
	{
		//We check if the value is less than the current
		//In essence checking for integer wrapping when 
		//the value tries to be negative in unsigned type
		while ((tempIndex - product) < tempIndex)
		{
			position.at(currentDimension) += 1;
			tempIndex -= product;
		}
		product /= m_piDimensions[currentDimension--];
	}
	return position;
}

UINT32 MatrixND::getIndexFromPosition(vector<UINT32> pos)
{
	if (isInMatrix(pos))
	{
		UINT32 index = 0;
		UINT32 product = 1;
		for (UINT16 i = 0; i < m_iDimensionality; i++)
		{
			product *= m_piDimensions[i];
		}
		for (UINT16 j = m_iDimensionality - 1; j >= 0; j--)
		{
			index += pos.at(j) * product;
			product /= pos.at(j);
		}
		return index;
	}
	else
	{
		//Returns the maximum UINT32 value as this is an extremely improbable scenario
		//This is only used because unsigned long int will not allow a positive value
		return UINT32_MAX;
	}
}

//-------------------------Operators------------------------------

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

MatrixND operator^(MatrixND mat1, MatrixND::OperatingDimensions_t dims)
{
	return MatrixND::transpose(mat1, dims);
}

//-----------------------------------------------------------------

void MatrixND::OperatingDimensions_t::set(UINT16 da, UINT16 db)
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
		this->db = db;
		return;
	}
}