#pragma once
#include<iostream>
#include <iomanip>
#include<functional>
#include<string>
#include<sstream>

namespace tnzr
{
	template <typename Tp_>
	class Matrix;

	template<typename Tp_>
	Matrix<Tp_> operator +(Matrix<Tp_> &first, Matrix<Tp_> &second);

	template<typename Tp_>
	Matrix<Tp_> operator -(Matrix<Tp_> &first, Matrix<Tp_> &second);

	template<typename Tp_>
	Matrix<Tp_> operator *(Matrix<Tp_> &first, Matrix<Tp_> &second);

	template<typename Tp_>
	Matrix<Tp_> operator *(const Matrix<Tp_> &first, const Tp_ scalar);

	template<typename Tp_>
	std::ostream& operator <<(std::ostream &out, Matrix<Tp_>& matrix);

	template <typename Tp_>
	class Matrix
	{
	protected:
		int id_, rows_, columns_;
		static int matrix_static_id_;
		Tp_* data_ptr_;

		struct RowWrapper 
		{
		private:
			Tp_* ptr_;
			int size_, id_;

		public:
			RowWrapper(Tp_ *ptr, int size, int id) : ptr_(ptr), size_(size), id_(id) {}

			Tp_ operator [](const int idx)  const
			{
				if (idx < size_ && idx >= 0)
				{
					return ptr_[idx];
				}
				else
				{
					throw MatrixException("Matrix::operator []: Wrong column index!", id_);
				}
			}

			Tp_& operator [](const int idx)
			{
				if (idx < size_ && idx >= 0)
				{
					return ptr_[idx];
				}
				else
				{
					throw MatrixException("Matrix::operator []: Wrong column index!", id_);
				}
			}
		};

	public: 

		class MatrixException : public std::exception
		{
		private:
			int id_info_ = -1;
			const char* msg_;

		public:

			MatrixException(const char* msg, int id) : id_info_(id), msg_(msg) {}

			MatrixException(const char* msg) :  msg_(msg) {}

			MatrixException() {}

			virtual const char* what() const
			{
				return msg_;
			}

			int which() const
			{
				return id_info_;
			}
		};

		Matrix(std::initializer_list<std::initializer_list<Tp_>> init_list);
		Matrix(int row, int col, std::function<Tp_()> customFunc = nullptr);
		Matrix(int row, int col, Tp_ data);
		Matrix(int size, std::function<Tp_()> func = nullptr);
		Matrix();
		~Matrix();
		Matrix(const Matrix &other);
		Matrix(Matrix &&other);

		Tp_ getValue(int idx_r, int idx_c) const;
		void setValue(int idx_r, int idx_c, Tp_ value);
		int getColNum() const;
		int getRowNum() const;
		bool checkMult(const Matrix &other) const;
		bool checkAdd(const Matrix &other) const;
		Tp_ max() const;
		Tp_ min() const;

		Matrix& operator =(const Matrix &other);
		Matrix& operator =(Matrix &&other);
		Matrix& operator +=(const Matrix &other);
		Matrix& operator -=(const Matrix &other);
		Matrix& operator *=(const Matrix &other);
		Matrix& operator *=(const Tp_ &scalar); 
		RowWrapper operator [](const int idx) const;
		friend Matrix operator + <>(Matrix &first,Matrix &second);
		friend Matrix operator - <>(Matrix &first, Matrix &second);
		friend Matrix operator * <>(Matrix &first, Matrix &second);
		friend Matrix operator * <>(const Matrix &first, const Tp_ scalar); 
		friend std::ostream& operator << <>(std::ostream &out, Matrix& matrix);

	};

	template <typename Tp_>
	int Matrix<Tp_>::matrix_static_id_ = 0;

	template<typename Tp_>
	Matrix<Tp_>::Matrix() : rows_(0), columns_(0), data_ptr_(nullptr), id_(++matrix_static_id_)
	{
		#ifdef _DEBUG
			std::cout << "Matrix " << this->id_ << " has been created" << std::endl;
		#endif
	}

	template<typename Tp_>
	Matrix<Tp_>::Matrix(std::initializer_list<std::initializer_list<Tp_>> init_list) 
	{
		columns_ = (int)(init_list.begin())->size();
		rows_ = (int)init_list.size();

		bool flag = true;
		for (int i = 1; i < rows_; i++)
			if ((init_list.begin() + i)->size() != columns_) flag = false;

		if (flag && columns_ != 0 && rows_ != 0)
		{

			try
			{
				this->data_ptr_ = new Tp_[rows_ * columns_];
			}
			catch (std::exception &ex)
			{
				delete[] data_ptr_;
				throw MatrixException("Matrix::Constructor: Bad memory allocation!");
			}

			for (int i = 0; i < rows_; i++)
				for (int j = 0; j < columns_; j++)
					data_ptr_[i * columns_ + j] = ((init_list.begin() + i)->begin())[j];

			id_ = ++this->static_id_;

			#ifdef _DEBUG
				std::cout << "Matrix " << this->id_ << " has been created" << std::endl;
			#endif
		}
		else
		{
			throw MatrixException("Matrix::Constructor: Bad initializer!");
		}
	}

	template<typename Tp_>
	Matrix<Tp_>::Matrix(int row, int col, std::function<Tp_()> customFunc) : rows_(row), columns_(col)
	{
		if (row > 0 && col > 0)
		{
			try
			{
				this->data_ptr_ = new Tp_[rows_ * columns_];
			}
			catch (std::exception &ex)
			{
				delete[] data_ptr_;
				throw MatrixException("Matrix::Constructor: Bad memory allocation!");
			}

			if (customFunc)
			{
				for (int i = 0; i < columns_*rows_; i++)
					data_ptr_[i] = customFunc();
			}

			id_ = ++matrix_static_id_;

			#ifdef _DEBUG
				std::cout << "Matrix " << this->id_ << " has been created" << std::endl;
			#endif
		}
		else
		{
			throw MatrixException("Matrix::Constructor: Wrong indices! Must be greater than zero!");
		}
	}

	template<typename Tp_>
	Matrix<Tp_>::Matrix(int row, int col,Tp_ data) : rows_(row), columns_(col)
	{
		if (row > 0 && col > 0)
		{
			try
			{
				this->data_ptr_ = new Tp_[rows_ * columns_];
			}
			catch (std::exception &ex)
			{
				delete[] data_ptr_;
				throw MatrixException("Matrix::Constructor: Bad memory allocation!");
			}

			for (int i = 0; i < columns_*rows_; i++)
				data_ptr_[i] = data;

			id_ = ++matrix_static_id_;

			#ifdef _DEBUG
				std::cout << "Matrix " << this->id_ << " has been created" << std::endl;
			#endif
		}
		else
		{
			throw MatrixException("Matrix::Constructor: Wrong indices! Must be greater than zero!");
		}
	}

	template<typename Tp_>
	Matrix<Tp_>::Matrix(int size, std::function<Tp_()> customFunc) : rows_(size), columns_(size)
	{
		if (size > 0)
		{
			try
			{
				this->data_ptr_ = new Tp_[rows_ * columns_];
			}
			catch (std::exception &ex)
			{
				delete[] data_ptr_;
				throw MatrixException("Matrix::Constructor: Bad memory allocation!");
			}

			if (customFunc)
			{
				for (int i = 0; i < columns_*rows_; i++)
					data_ptr_[i] = customFunc();
			}

			this->id_ = ++this->static_id_;

			#ifdef _DEBUG
				std::cout << "Matrix " << this->id_ << " has been created" << std::endl;
			#endif
		}
		else
		{
			throw MatrixException("Matrix::Constructor: Wrong indices! Must be greater than zero!");
		}
	}

	template<typename Tp_>
	Matrix<Tp_>::Matrix(const Matrix &other) : columns_(other.columns_), rows_(other.rows_)
	{
		try
		{
			this->data_ptr_ = new Tp_[rows_ * columns_];
		}
		catch (std::exception &ex)
		{
			delete[] data_ptr_;
			throw MatrixException("Matrix::CopyConstructor: Bad memory allocation!");
		}

		for (int i = 0; i < columns_ * rows_; i++)
			this->data_ptr_[i] = other.data_ptr_[i];

		id_ = ++matrix_static_id_;

		#ifdef _DEBUG
			std::cout << "Matrix " << this->id_ << " has been created" << std::endl;
		#endif
	}

	template<typename Tp_>
	Matrix<Tp_>::Matrix(Matrix &&other) : columns_(other.columns_), rows_(other.rows_), data_ptr_(other.data_ptr_)
	{
		other.columns_ = 0;
		other.rows_ = 0;
		other.id_ = 0;
		other.data_ptr_ = nullptr;

		#ifdef _DEBUG
			std::cout << "Move assigned" << std::endl;
		#endif
	}

	template<typename Tp_>
	Matrix<Tp_>& Matrix<Tp_>::operator =(Matrix<Tp_> &&other)
	{
		if (this != &other)
		{
			delete[] this->data_ptr_;
			this->columns_ = other.columns_;
			this->rows_ = other.rows_;
			this->data_ptr_ = other.data_ptr_;
			other.columns_ = 0;
			other.rows_ = 0;
			other.id_ = 0;
			other.data_ptr_ = nullptr;

			#ifdef _DEBUG
				std::cout << "Move assigned" << std::endl;
			#endif

			return *this;
		}
		else
		{
			throw MatrixException("Matrix::move operator: Can't be assigned to itself!", id_);
		}

	}

	template<typename Tp_>
	Matrix<Tp_>::~Matrix()
	{
		delete[] data_ptr_;

		#ifdef _DEBUG
			std::cout << "Matrix " << this->id_ << " has been destroyed" << std::endl;
		#endif
	}

	template<typename Tp_>
	Tp_ Matrix<Tp_>::getValue(int idx_r, int idx_c) const
	{
		if (idx_r >= 0 && idx_c >= 0 && idx_r < rows_ && idx_c < columns_)
		{
			return data_ptr_[idx_r * columns_ + idx_c];
		}
		else
		{
			throw MatrixException("Matrix::getValue: Wrong indices! Must be not less than zero and less than size of matrix!", id_);
		}
	}

	template<typename Tp_>
	void Matrix<Tp_>::setValue(int idx_r, int idx_c, Tp_ value)
	{
		if (idx_r >= 0 && idx_c >= 0 && idx_r < rows_ && idx_c < columns_)
		{
			data_ptr_[idx_r * columns_ + idx_c] = value;
		}
		else
		{
			throw MatrixException("Matrix::setValue: Wrong indices! Must be not less than zero and less than size of matrix!", id_);
		}
	}

	template<typename Tp_>
	int Matrix<Tp_>::getColNum() const { return columns_; }

	template<typename Tp_>
	int Matrix<Tp_>::getRowNum() const { return rows_; }

	template<typename Tp_>
	bool Matrix<Tp_>::checkMult(const Matrix &other) const { return this->columns_ == other.rows_; }

	template<typename Tp_>
	bool Matrix<Tp_>::checkAdd(const Matrix &other) const { return this->columns_ == other.columns_ && this->rows_ == other.rows_; }

	template<typename Tp_>
	Tp_ Matrix<Tp_>::max() const
	{
		if (data_ptr_ != nullptr)
		{
			Tp_ max = data_ptr_[0];

			for (int i = 0; i < columns_ * rows_; i++)
				if (data_ptr_[i] > max) max = data_ptr_[i];

			return max;
		}
		else
		{
			throw MatrixException("Matrix::max: Matrix is invalid!" , id_);
		}
	}

	template<typename Tp_>
	Tp_ Matrix<Tp_>::min() const
	{
		if (data_ptr_ != nullptr)
		{
			Tp_ min = data_ptr_[0];

			for (int i = 0; i < columns_ * rows_; i++)
				if (data_ptr_[i] < min) min = data_ptr_[i];

			return min;
		}
		else
		{
			throw MatrixException("Matrix::min: Matrix is invalid!", id_);
		}
	}

	template<typename Tp_>
	Matrix<Tp_>& Matrix<Tp_>::operator =(const Matrix<Tp_> &other)
	{
		delete[] this->data_ptr_;
		this->columns_ = other.columns_;
		this->rows_ = other.rows_;

		try
		{
			this->data_ptr_ = new Tp_[rows_ * columns_];
		}
		catch (std::exception &ex)
		{
			delete[] data_ptr_;
			throw MatrixException("Matrix::operator = : Bad memory allocation!", id_);
		}

		for (int i = 0; i < columns_ * rows_; i++)
			this->data_ptr_[i] = other.data_ptr_[i];

		return *this;
	}

	template<typename Tp_>
	Matrix<Tp_>& Matrix<Tp_>::operator +=(const Matrix<Tp_> &other)
	{
		if (this->checkAdd(other))
		{
			for (int i = 0; i < columns_ * rows_; i++)
				this->data_ptr_[i] += other.data_ptr_[i];

			return *this;
		} 
		else
		{
			throw MatrixException("Matrix::operator +=: Matrix can't be added!", other.id_);
		}
	}

	template<typename Tp_>
	Matrix<Tp_>& Matrix<Tp_>::operator -=(const Matrix<Tp_> &other)
	{
		if (this->checkAdd(other))
		{
			for (int i = 0; i < columns_ * rows_; i++)
				this->data_ptr_[i] -= other.data_ptr_[i];

			return *this;
		}
		else
		{
			throw MatrixException("Matrix::operator -=: Matrix can't be substracted!", other.id_);
		}
	}

	template<typename Tp_>
	Matrix<Tp_>& Matrix<Tp_>::operator *=(const Matrix<Tp_> &other)
	{
		if (this->checkMult(other))
		{
			Tp_ *tmp = nullptr;

			try
			{
				tmp = new Tp_[this->rows_ * other.columns_];
			}
			catch (std::exception &ex)
			{
				delete[] data_ptr_;
				throw MatrixException("Matrix::operator *= : Bad memory allocation!", id_);
			}
						
			try
			{
				for (int i = 0; i < this->rows_; i++)
				{
					for (int j = 0; j < other.columns_; j++)
					{
						tmp[i * other.columns_ + j] = this->data_ptr_[i * this->columns_] * other.data_ptr_[j];
						for (int k = 1; k < this->columns_; k++)
						{
							tmp[i * other.columns_ + j] += this->data_ptr_[i * this->columns_ + k] * other.data_ptr_[k * other.columns_ + j];
						}
					}
				}
				
			}
			catch (std::exception &ex)
			{
				delete[] tmp;
				throw ex;
			}

			delete[] this->data_ptr_;
			this->columns_ = other.columns_;
			this->data_ptr_ = tmp;

			return *this;
		}
		else
		{
			throw MatrixException("Matrix::operator *=: Matrix can't be multiplier!", other.id_);
		}
	}

	template<typename Tp_>
	Matrix<Tp_>& Matrix<Tp_>::operator *=(const Tp_ &scalar)
	{
		for (int i = 0; i < columns_ * rows_; i++)
			this->data_ptr_[i] *= scalar;

		return *this;
	}

	template<typename Tp_>
	typename Matrix<Tp_>::RowWrapper Matrix<Tp_>::operator [](const int idx) const
	{
		if (idx < rows_ && idx >= 0)
		{
			RowWrapper new_wrapper(&(data_ptr_[idx * columns_]), columns_, id_);
			return new_wrapper;
		}
		else
		{
			throw MatrixException("Matrix::operator []: Wrong row index!", id_);
		}
	}

	template<typename Tp_>
	Matrix<Tp_> operator +(const Matrix<Tp_> &first, const Matrix<Tp_> &second) 
	{
		if (first.checkAdd(second))
		{
			try
			{
				Matrix<Tp_> new_matrix(first.rows_, first.columns_);

				for (int i = 0; i < new_matrix.columns_ * new_matrix.rows_; i++)
					new_matrix.data_ptr_[i] = first.data_ptr_[i] + second.data_ptr_[i];

				return new_matrix;
			}
			catch (typename Matrix<Tp_>::MatrixException &ex)
			{
				throw ex;
			}
		}
		else
		{
			throw typename Matrix<Tp_>::MatrixException("Matrix::operator +: Matrix can't be added!", second.id_);
		}
	}

	template<typename Tp_>
	Matrix<Tp_> operator -(Matrix<Tp_> &first, Matrix<Tp_> &second)
	{
		if (first.checkAdd(second))
		{
			try
			{
				Matrix<Tp_> new_matrix(first.rows_, first.columns_);

				for (int i = 0; i < new_matrix.columns_ * new_matrix.rows_; i++)
					new_matrix.data_ptr_[i] = first.data_ptr_[i] - second.data_ptr_[i];

				return new_matrix;
			}
			catch (typename Matrix<Tp_>::MatrixException &ex)
			{
				throw ex;
			}
		}
		else
		{
			throw typename Matrix<Tp_>::MatrixException("Matrix::operator -: Matrix can't be substracted!", second.id_);
		}
	}

	template<typename Tp_>
	Matrix<Tp_> operator *(Matrix<Tp_> &first, Matrix<Tp_> &second)
	{
		if (first.checkMult(second))
		{
			try
			{
				Matrix<Tp_> new_matrix(first.rows_, second.columns_);

				for (int i = 0; i < first.rows_; i++)
				{
					for (int j = 0; j < second.columns_; j++)
					{
						new_matrix.data_ptr_[i * second.columns_ + j] = first.data_ptr_[i * first.columns_] * second.data_ptr_[j];
						for (int k = 1; k < first.columns_; k++)
						{
							new_matrix.data_ptr_[i * second.columns_ + j] += first.data_ptr_[i * first.columns_ + k] *
								second.data_ptr_[k * second.columns_ + j];
						}
					}
				}

				return new_matrix;
			}
			catch (typename Matrix<Tp_>::MatrixException &ex)
			{
				throw ex;
			}
		}
		else
		{
			throw typename Matrix<Tp_>::MatrixException("Matrix::operator *=: Matrix can't be multiplier!", second.id_);
		}
	}

	template<typename Tp_>
	Matrix<Tp_> operator *(const Matrix<Tp_> &first, const Tp_ scalar)
	{
		try
		{
			Matrix<Tp_> new_matrix(first.rows_, first.columns_);

			for (int i = 0; i < first.columns_ * first.rows_; i++)
				new_matrix.data_ptr_[i] = first.data_ptr_[i] * scalar;

			return new_matrix;
		}
		catch (typename Matrix<Tp_>::MatrixException &ex)
		{
			throw ex;
		}
	}

	template<typename Tp_>
	std::ostream& operator <<(std::ostream &out, Matrix<Tp_>& matrix) 
	{
		int width = out.width();
		int precision = out.precision();

		for (int i = 0; i < matrix.rows_; i++)
		{
			for (int j = 0; j < matrix.columns_; j++)
			{
				out << std::setfill(' ') << std::setw(width)  << std::setprecision(precision) << matrix.data_ptr_[i * matrix.columns_ + j] << " ";
			}
			out << std::endl;
		}
		
		return out;
	}

}

