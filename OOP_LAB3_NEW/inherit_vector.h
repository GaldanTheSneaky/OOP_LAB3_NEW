#pragma once
#include"matrix.h"
#define _USE_MATH_DEFINES
#include<math.h>
#include<iostream>

namespace tnzr
{
	template <typename Tp_>
	class Vector;

	template<typename Tp_>
	Vector<Tp_> crossProduct(const Vector<Tp_>& a, const Vector<Tp_>& b);

	template<typename Tp_>
	Tp_ dotProduct(const Vector<Tp_>& a, const Vector<Tp_>& b);

	template<typename Tp_>
	Tp_ cos(const Vector<Tp_>& a, const Vector<Tp_>& b);

	template<typename Tp_>
	Tp_ sin(const Vector<Tp_>& a, const Vector<Tp_>& b);

	template<typename Tp_>
	Tp_ getAngle(Vector<Tp_>& a, Vector<Tp_>& b);

	//template<typename Tp_>
	//Vector<Tp_> operator +(Vector<Tp_>& first, Vector<Tp_>& second);

	//template<typename Tp_>
	//Vector<Tp_> operator -(Vector<Tp_>& first, Vector<Tp_>& second);

	//template<typename Tp_>
	//Vector<Tp_> operator *(Vector<Tp_>& first, Vector<Tp_>& second);

	//template<typename Tp_>
	//Vector<Tp_> operator *(const Vector<Tp_>& first, const Tp_ scalar);

	//template<typename Tp_>
	//std::ostream& operator <<(std::ostream& out, Vector<Tp_>& matrix);

	template<typename Tp_>
	class Vector : public Matrix < Tp_ >
	{
	private:
		static int vector_static_id_;
	public:
		Vector(std::initializer_list<Tp_> init_list);
		Vector(int size, std::function<Tp_()> customFunc = nullptr);
		Vector(int size, Tp_ data);
		Vector();

		Tp_ getMagnitude() const;
		Vector<Tp_> getScaled(const double scalar) const;
		Vector<Tp_> getNormalized() const;
		void scale(const Tp_& scalar);
		void normalize();


		template<typename Tp_>
		friend Vector<Tp_> crossProduct(const Vector<Tp_>& a, const Vector<Tp_>& b);

		template<typename Tp_>
		friend Tp_ dotProduct(const Vector<Tp_>& a, const Vector<Tp_>& b);

		template<typename Tp_>
		friend Tp_ cos(const Vector<Tp_>& a, const Vector<Tp_>& b);

		template<typename Tp_>
		friend Tp_ sin(const Vector<Tp_>& a, const Vector<Tp_>& b);

		template<typename Tp_>
		friend Tp_ getAngle(Vector<Tp_>& a, Vector<Tp_>& b);

		/*Vector<Tp_>& operator =(const Vector<Tp_>& other);*/

		//template<typename Tp_>
		//friend Vector<Tp_> operator +(const Vector<Tp_>& first, const Vector<Tp_>& second)
		//{
		//	return static_cast<Matrix<Tp_>&>(first) + static_cast<Matrix<Tp_>&>(second);
		//	//return operator +(first, second);
		//}
		//friend Matrix operator - <>(Matrix& first, Matrix& second);
		//friend Matrix operator * <>(Matrix& first, Matrix& second);
		//friend Matrix operator * <>(const Matrix& first, const Tp_ scalar);

		Tp_ operator [](const int idx) const;
		Tp_& operator [](const int idx);


	};

	template <typename Tp_>
	int Vector<Tp_>::vector_static_id_ = 0;

	template <typename Tp_>
	Vector<Tp_>::Vector(std::initializer_list<Tp_> init_list)
	{
		this->columns_ = 1;
		this->rows_ = (int)init_list.size();

		if (this->rows_ != 0)
		{
			try
			{
				this->data_ptr_ = new Tp_[this->rows_ * this->columns_];
			}
			catch (std::exception &ex)
			{
				delete[] this->data_ptr_;
				throw Matrix<Tp_>::MatrixException("Matrix::Constructor: Bad memory allocation!");
			}

			for (int i = 0; i < this->rows_; i++)
				this->data_ptr_[i] = init_list.begin()[i];

			this->id_ = ++vector_static_id_;

#ifdef _DEBUG
			std::cout << "Matrix " << this->id_ << " has been created" << std::endl;
#endif
		}
		else
		{
			throw Matrix<Tp_>::MatrixException("Matrix::Constructor: Bad initializer!");
		}
	}

	template<typename Tp_>
	Vector<Tp_>::Vector() : Matrix<Tp_>() {}

	template<typename Tp_>
	Vector<Tp_>::Vector(int size, std::function<Tp_()> customFunc) : Matrix<Tp_>(size, 1, customFunc){}

	template<typename Tp_>
	Vector<Tp_>::Vector(int size, Tp_ data) : Matrix<Tp_>(size, 1, data){}

	template<typename Tp_>
	Tp_ Vector<Tp_>::getMagnitude() const
	{
		Tp_ sum = 0;
		for (int i = 0; i < this->rows_; i++)
			sum += this->data_ptr_[i] * this->data_ptr_[i];

		return sqrt(sum);
	}

	template<typename Tp_>
	void Vector<Tp_>::scale(const Tp_& scalar)
	{
		for (int i = 0; i < this->rows_; i++)
			this->data_ptr_[i] *= scalar;
	}

	template<typename Tp_>
	Vector<Tp_> Vector<Tp_>::getScaled(const double scalar) const
	{
		Vector newVec = *this;
		newVec.scale(scalar);

		return newVec;
	}

	template<typename Tp_>
	void Vector<Tp_> ::normalize()
	{
		Tp_ module = this->getMagnitude();
		if (module != 0)
		{
			for (int i = 0; i < this->rows_; i++)
				this->data_ptr_[i] /= module;
		}
	}

	template<typename Tp_>
	Vector<Tp_> Vector<Tp_>::getNormalized() const
	{
		Vector newVec = *this;
		newVec.normalize();

		return newVec;
	}

	template<typename Tp_>
	Tp_ Vector<Tp_>::operator [](const int idx)  const
	{
		if (idx < this->rows_ && idx >= 0)
		{
			return this->data_ptr_[idx];
		}
		else
		{
			throw MatrixException("Matrix::operator []: Wrong column index!", this->id_);
		}
	}

	template<typename Tp_>
	Tp_& Vector<Tp_>::operator [](const int idx)
	{
		if (idx < this->rows_ && idx >= 0)
		{
			return this->data_ptr_[idx];
		}
		else
		{
			throw MatrixException("Matrix::operator []: Wrong column index!", this->id_);
		}
	}

	template<typename Tp_>
	Vector<Tp_> crossProduct(const Vector<Tp_>& a, const Vector<Tp_>& b)
	{
		if (a.rows_ == 3 && b.rows_ == 3)
		{
			Vector<Tp_> newVec(3);
			newVec[0] = a[1] * b[2] - a[2] * b[1];
			newVec[1] = a[2] * b[0] - a[0] * b[2];
			newVec[2] = a[0] * b[1] - a[1] * b[0];

			return newVec;
		}
		throw Matrix<Tp_>::MatrixException("Matrix::Must be 3 deminsional");
	}

	template<typename Tp_>
	Tp_ dotProduct(const Vector<Tp_>& a, const Vector<Tp_>& b)
	{
		Tp_ sum = 0;
		if (a.rows_ == b.rows_)
		{
			for (int i = 0; i < a.rows_; i++)
				sum += a[i] * b[i];

			return sum;
		}
		throw Matrix<Tp_>::MatrixException("Matrix::Vectors's dimensions must be equal");
	}

	template<typename Tp_>
	Tp_ cos(const Vector<Tp_>& a, const Vector<Tp_>& b)
	{
		return dotProduct(a, b) / (a.getMagnitude() * b.getMagnitude());

	}

	template<typename Tp_>
	Tp_ sin(const Vector<Tp_>& a, const Vector<Tp_>& b)
	{
		return crossProduct(a, b).getMagnitude() / (a.getMagnitude() * b.getMagnitude());
	}

	template<typename Tp_>
	Tp_ getAngle(Vector<Tp_>& a, Vector<Tp_>& b)
	{
		return atan2(sin(a, b), cos(a, b)) * 180 / M_PI;
	}

	//template<typename Tp_>
	//Vector<Tp_>& Vector<Tp_>::operator =(const Vector<Tp_>& other)
	//{
	//	Matrix<Tp_>::operator=(other);
	//	return *this;
	//}
	//Vector& Vector<Tp_>operator =(Vector&& other);

	//template<typename Tp_>
	//Vector<Tp_> operator +(const Vector<Tp_>& first, const Vector<Tp_>& second)
	//{
	//	return static_cast<Matrix<Tp_>&>(first) + static_cast<Matrix<Tp_>&>(second);
	//}

	//Vector& Vector<Tp_>operator -=(const Vector& other);
	//Vector& Vector<Tp_>operator *=(const Vector& other);
	//Vector& Vector<Tp_>operator *=(const Tp_& scalar);
}


