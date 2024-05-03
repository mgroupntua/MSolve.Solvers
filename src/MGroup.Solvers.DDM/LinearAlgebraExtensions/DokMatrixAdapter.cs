namespace MGroup.Solvers.DDM.LinearAlgebraExtensions
{
	using System;

	using MGroup.LinearAlgebra.Matrices;
	using MGroup.LinearAlgebra.Providers;
	using MGroup.LinearAlgebra.Reduction;
	using MGroup.LinearAlgebra.Vectors;

	public class DokMatrixAdapter<TMatrix> : IMatrix where TMatrix : IIndexable2D
	{
		public DokMatrixAdapter(TMatrix dok)
		{
			this.DOK = dok;
		}

		public TMatrix DOK { get; set; }

		public int NumColumns => DOK.NumColumns;

		public int NumRows => DOK.NumRows;

		public MatrixSymmetry MatrixSymmetry => MatrixSymmetry.Unknown;

		public double this[int rowIdx, int colIdx] => DOK[rowIdx, colIdx];

		public IMatrix Axpy(IMatrixView otherMatrix, double otherCoefficient)
		{
			throw new NotImplementedException();
		}

		public void AxpyIntoThis(IMatrixView otherMatrix, double otherCoefficient)
		{
			throw new NotImplementedException();
		}

		public void Clear() => throw new NotImplementedException();

		public IMatrix Copy(bool copyIndexingData = false)
		{
			throw new NotImplementedException();
		}

		public Matrix CopyToFullMatrix()
		{
			var full = Matrix.CreateZero(NumRows, NumColumns);
			for (int j = 0; j < full.NumColumns; ++j)
			{
				for (int i = 0; i < full.NumRows; ++i)
				{
					full[i, j] = DOK[i, j];
				}
			}
			return full;
		}

		public IMatrix DoEntrywise(IMatrixView matrix, Func<double, double, double> binaryOperation)
		{
			throw new NotImplementedException();
		}

		public void DoEntrywiseIntoThis(IMatrixView matrix, Func<double, double, double> binaryOperation)
		{
			throw new NotImplementedException();
		}

		public IMatrix DoToAllEntries(Func<double, double> unaryOperation)
		{
			throw new NotImplementedException();
		}

		public void DoToAllEntriesIntoThis(Func<double, double> unaryOperation)
		{
			throw new NotImplementedException();
		}

		public bool Equals(IIndexable2D other, double tolerance = 1E-13) => DOK.Equals(other, tolerance);

		public Vector GetColumn(int colIndex) => throw new NotImplementedException();

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

		public void LinearCombinationIntoThis(double thisCoefficient, IMatrixView otherMatrix, double otherCoefficient)
		{
			throw new NotImplementedException();
		}

		public IVector Multiply(IVectorView vector, bool transposeThis = false)
		{
			throw new NotImplementedException();
		}

		public void MultiplyIntoResult(IVectorView lhsVector, IVector rhsVector, bool transposeThis = false)
		{
			throw new NotImplementedException();
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

		public void ScaleIntoThis(double scalar)
		{
			throw new NotImplementedException();
		}

		public void SetEntryRespectingPattern(int rowIdx, int colIdx, double value)
		{
			throw new NotImplementedException();
		}

		public IMatrix Transpose()
		{
			throw new NotImplementedException();
		}
	}
}
