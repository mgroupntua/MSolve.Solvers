namespace MGroup.Solvers.DDM.LinearAlgebraExtensions
{
	using System;

	using MGroup.LinearAlgebra.Commons;
	using MGroup.LinearAlgebra.Matrices;
	using MGroup.LinearAlgebra.Providers;
	using MGroup.LinearAlgebra.Reduction;
	using MGroup.LinearAlgebra.Vectors;

	public class CsrMatrixFaster : IMatrixView
	{
		private readonly CsrMatrix simpleCsr;
		private readonly double[] values;
		private readonly int[] colIndices;
		private readonly int[] rowOffsets;

		public CsrMatrixFaster(CsrMatrix matrix)
		{
			this.simpleCsr = matrix;
			this.NumRows = matrix.NumRows;
			this.NumColumns = matrix.NumColumns;
			this.values = matrix.RawValues;
			this.colIndices = matrix.RawColIndices;
			this.rowOffsets = matrix.RawRowOffsets;
		}

		public double this[int rowIdx, int colIdx] => throw new NotImplementedException();

		public int NumColumns { get; }

		public int NumRows { get; }

		public double[] RawValues => values;
		public int[] RawColIndices => colIndices;
		public int[] RawRowOffsets => rowOffsets;

		public CsrMatrix AsSimpleCsr => simpleCsr;

		public MatrixSymmetry MatrixSymmetry => MatrixSymmetry.Unknown;

		public static Vector operator *(CsrMatrixFaster matrixLeft, Vector vectorRight)
			=> matrixLeft.Multiply(vectorRight, false);

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
			var result = Vector.CreateZero(transposeThis ? NumColumns : NumRows);
			MultiplyIntoResult(vector, result, transposeThis);
			return result;
		}

		public Vector Multiply(Vector vector, bool transposeThis = false)
		{
			var result = Vector.CreateZero(transposeThis ? NumColumns : NumRows);
			MultiplyIntoResult(vector, result, transposeThis);
			return result;
		}

		public void MultiplyIntoResult(IVectorView lhsVector, IVector rhsVector, bool transposeThis = false)
		{
			MultiplyIntoResult((Vector)lhsVector, (Vector)rhsVector, transposeThis);
		}

		public void MultiplyIntoResult(Vector lhsVector, Vector rhsVector, bool transposeThis = false)
		{
			double[] x = lhsVector.RawData;
			double[] y = rhsVector.RawData;
			if (transposeThis)
			{
				Preconditions.CheckMultiplicationDimensions(NumRows, x.Length);
				Preconditions.CheckSystemSolutionDimensions(NumColumns, y.Length);
				CsrTransposedTimesVectorUnsafe(NumRows, NumColumns, values, rowOffsets, colIndices, x, y);
			}
			else
			{
				Preconditions.CheckMultiplicationDimensions(NumColumns, x.Length);
				Preconditions.CheckSystemSolutionDimensions(NumRows, y.Length);
				CsrTimesVectorUnsafe(NumRows, NumColumns, values, rowOffsets, colIndices, x, y);
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

		public IMatrix Scale(double scalar)
		{
			throw new NotImplementedException();
		}

		public IMatrix Transpose()
		{
			throw new NotImplementedException();
		}

		private static void CsrTimesVector(int numRowsA, int numColsA, double[] valuesA, int[] rowOffsetsA,
			int[] colIndicesA, double[] x, int offsetX, double[] y, int offsetY)
		{
			Array.Clear(y, offsetY, numRowsA);
			for (int i = 0; i < numRowsA; ++i)
			{
				double dot = 0.0;
				int rowStart = rowOffsetsA[i]; //inclusive
				int rowEnd = rowOffsetsA[i + 1]; //exclusive
				for (int k = rowStart; k < rowEnd; ++k) dot += valuesA[k] * x[offsetX + colIndicesA[k]];
				y[offsetY + i] = dot;
			}
		}

		private static void CsrTimesVectorNoOffset(int numRowsA, int numColsA, double[] valuesA, int[] rowOffsetsA,
			int[] colIndicesA, double[] x, double[] y)
		{
			//Array.Clear(y, 0, numRowsA);
			for (int i = 0; i < numRowsA; ++i)
			{
				double dot = 0.0;
				int rowStart = rowOffsetsA[i]; //inclusive
				int rowEnd = rowOffsetsA[i + 1]; //exclusive
				for (int k = rowStart; k < rowEnd; ++k) dot += valuesA[k] * x[colIndicesA[k]];
				y[i] = dot;
			}
		}

		private static void CsrTimesVectorUnsafe(int numRowsA, int numColsA, double[] valuesA, int[] rowOffsetsA,
			int[] colIndicesA, double[] x, double[] y)
		{
			unsafe
			{
				fixed (double* values = &valuesA[0])
				fixed (int* rowOffsets = &rowOffsetsA[0])
				fixed (int* colIndices = &colIndicesA[0])
				fixed (double* lhs = &x[0])
				fixed (double* rhs = &y[0])
				{
					for (int i = 0; i < numRowsA; ++i)
					{
						double dot = 0.0;
						int rowStart = rowOffsets[i]; //inclusive
						int rowEnd = rowOffsets[i + 1]; //exclusive
						for (int k = rowStart; k < rowEnd; ++k) dot += values[k] * lhs[colIndicesA[k]];
						rhs[i] = dot;
					}

				}
			}
		}

		private static void CsrTransposedTimesVector(int numRowsA, int numColsA, double[] valuesA, int[] rowOffsetsA,
			int[] colIndicesA, double[] x, int offsetX, double[] y, int offsetY)
		{
			Array.Clear(y, offsetY, numColsA);
			// A^T * x = linear combination of columns of A^T = rows of A, with the entries of x as coefficients
			for (int i = 0; i < numRowsA; ++i)
			{
				double scalar = x[offsetX + i];
				int rowStart = rowOffsetsA[i]; //inclusive
				int rowEnd = rowOffsetsA[i + 1]; //exclusive
				for (int k = rowStart; k < rowEnd; ++k)
				{
					y[offsetY + colIndicesA[k]] += scalar * valuesA[k];
				}
			}
		}

		private static void CsrTransposedTimesVectorNoOffset(int numRowsA, int numColsA, double[] valuesA, int[] rowOffsetsA,
			int[] colIndicesA, double[] x, double[] y)
		{
			Array.Clear(y, 0, numColsA);
			// A^T * x = linear combination of columns of A^T = rows of A, with the entries of x as coefficients
			for (int i = 0; i < numRowsA; ++i)
			{
				double scalar = x[i];
				int rowStart = rowOffsetsA[i]; //inclusive
				int rowEnd = rowOffsetsA[i + 1]; //exclusive
				for (int k = rowStart; k < rowEnd; ++k)
				{
					y[colIndicesA[k]] += scalar * valuesA[k];
				}
			}
		}

		private static void CsrTransposedTimesVectorUnsafe(int numRowsA, int numColsA, double[] valuesA, int[] rowOffsetsA,
			int[] colIndicesA, double[] x, double[] y)
		{
			Array.Clear(y, 0, numColsA);
			unsafe
			{
				fixed (double* values = &valuesA[0])
				fixed (int* rowOffsets = &rowOffsetsA[0])
				fixed (int* colIndices = &colIndicesA[0])
				fixed (double* lhs = &x[0])
				fixed (double* rhs = &y[0])
				{
					// A^T * x = linear combination of columns of A^T = rows of A, with the entries of x as coefficients
					for (int i = 0; i < numRowsA; ++i)
					{
						double scalar = lhs[i];
						int rowStart = rowOffsetsA[i]; //inclusive
						int rowEnd = rowOffsetsA[i + 1]; //exclusive
						for (int k = rowStart; k < rowEnd; ++k)
						{
							rhs[colIndicesA[k]] += scalar * valuesA[k];
						}
					}
				}
			}
		}
	}
}
