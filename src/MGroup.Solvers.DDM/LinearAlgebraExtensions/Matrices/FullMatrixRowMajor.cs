namespace MGroup.Solvers.DDM.LinearAlgebraExtensions.Matrices
{
	using System;
	using System.Diagnostics;

	using MGroup.LinearAlgebra.Matrices;
	using MGroup.LinearAlgebra.Providers;
	using MGroup.LinearAlgebra.Reduction;
	using MGroup.LinearAlgebra.Vectors;

	public class FullMatrixRowMajor : IMatrixView
	{
		private readonly double[] data;

		private FullMatrixRowMajor(int numRows, int numColumns, double[] data)
		{
			this.NumRows = numRows;
			this.NumColumns = numColumns;
			this.data = data;
		}

		public double this[int rowIdx, int colIdx] => throw new NotImplementedException();

		public int NumColumns { get; }

		public int NumRows { get; }

		public double[] RawData => data;

		public MatrixSymmetry MatrixSymmetry => MatrixSymmetry.Unknown;

		public static FullMatrixRowMajor CreateFromArray(int numRows, int numColumns, double[] data)
		{
			return new FullMatrixRowMajor(numRows, numColumns, data);
		}

		public static FullMatrixRowMajor CreateFromZero(int numRows, int numColumns)
		{
			return new FullMatrixRowMajor(numRows, numColumns, new double[numRows * numColumns]);
		}

		public IMatrix Axpy(IMatrixView otherMatrix, double otherCoefficient)
		{
			throw new NotImplementedException();
		}

		public IMatrix Copy(bool copyIndexingData = false)
		{
			throw new NotImplementedException();
		}

		public Matrix CopyToFullMatrix()
		{
			throw new NotImplementedException();
		}

		public IMatrix DoEntrywise(IMatrixView matrix, Func<double, double, double> binaryOperation)
		{
			throw new NotImplementedException();
		}

		public IMatrix DoToAllEntries(Func<double, double> unaryOperation)
		{
			throw new NotImplementedException();
		}

		public bool Equals(IIndexable2D other, double tolerance = 1E-13)
		{
			throw new NotImplementedException();
		}

		public Vector GetColumn(int colIndex)
		{
			throw new NotImplementedException();
		}

		public Vector GetRow(int rowIndex)
		{
			throw new NotImplementedException();
		}

		public IMatrix GetSubmatrix(int[] rowIndices, int[] colIndices)
		{
			throw new NotImplementedException();
		}

		public IMatrix GetSubmatrix(int rowStartInclusive, int rowEndExclusive, int colStartInclusive, int colEndExclusive)
		{
			throw new NotImplementedException();
		}

		public IMatrix LinearCombination(double thisCoefficient, IMatrixView otherMatrix, double otherCoefficient)
		{
			throw new NotImplementedException();
		}

		public IVector Multiply(IVectorView vector, bool transposeThis = false)
		{
			if (vector is Vector denseVector)
			{
				return Multiply(denseVector, transposeThis);
			}
			else
			{
				throw new NotImplementedException();
			}
		}

		public Vector Multiply(Vector vector, bool transposeThis = false)
		{
			if (transposeThis)
			{
				Vector result = Vector.CreateZero(this.NumColumns);
				TransposeMultiplyVector(vector, result);
				return result;
			}
			else
			{
				Vector result = Vector.CreateZero(this.NumRows);
				MultiplyVector(vector, result);
				return result;
			}
		}

		public void MultiplyIntoResult(IVectorView lhsVector, IVector rhsVector, bool transposeThis = false)
		{
			if ((lhsVector is Vector denseLhs) && (rhsVector is Vector denseRhs))
			{
				MultiplyIntoResult(denseLhs, denseRhs, transposeThis);
			}
			else
			{
				throw new NotImplementedException();
			}
		}

		public void MultiplyIntoResult(Vector lhsVector, Vector rhsVector, bool transposeThis = false)
		{
			if (transposeThis)
			{
				TransposeMultiplyVector(lhsVector, rhsVector);
			}
			else
			{
				MultiplyVector(lhsVector, rhsVector);
			}
		}

		public Matrix MultiplyLeft(IMatrixView other, bool transposeThis = false, bool transposeOther = false)
		{
			throw new NotImplementedException();
		}

		public Matrix MultiplyRight(IMatrixView other, bool transposeThis = false, bool transposeOther = false)
		{
			throw new NotImplementedException();
		}

		public double Reduce(double identityValue, ProcessEntry processEntry, ProcessZeros processZeros, Finalize finalize)
		{
			throw new NotImplementedException();
		}

		public void SetRow(int rowIdx, Vector rowValues)
		{
			Debug.Assert(rowIdx >= 0 && rowIdx < this.NumRows);
			Debug.Assert(rowValues.Length == this.NumColumns);
			Array.Copy(rowValues.RawData, 0, this.data, rowIdx * NumColumns, NumColumns);
		}

		public IMatrix Scale(double scalar)
		{
			throw new NotImplementedException();
		}

		public IMatrix Transpose()
		{
			throw new NotImplementedException();
		}

		private void MultiplyVector(Vector lhs, Vector rhs)
		{
			Debug.Assert(lhs.Length == this.NumColumns);
			Debug.Assert(rhs.Length == this.NumRows);
			double[] x = lhs.RawData;
			double[] y = rhs.RawData;
			for (int i = 0; i < NumRows; ++i)
			{
				int offset = i * NumColumns;
				double sum = 0;
				for (int j = 0; j < NumColumns; ++j)
				{
					sum += data[offset + j] * x[j];
				}
				y[i] = sum;
			}
		}

		private void TransposeMultiplyVector(Vector lhs, Vector rhs)
		{

			Debug.Assert(lhs.Length == this.NumRows);
			Debug.Assert(rhs.Length == this.NumColumns);
			double[] x = lhs.RawData;
			double[] y = rhs.RawData;
			for (int i = 0; i < NumColumns; ++i)
			{
				double sum = 0;
				for (int j = 0; j < NumRows; ++j)
				{
					sum += data[i * NumRows + j] * x[j];

				}
				y[i] = sum;
			}
		}
	}
}
