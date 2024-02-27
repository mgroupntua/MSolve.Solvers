using System;
using System.Collections.Generic;
using System.Text;
using MGroup.LinearAlgebra.Commons;
using MGroup.LinearAlgebra.Matrices;
using MGroup.LinearAlgebra.Vectors;

namespace MGroup.Solvers.DDM.Mappings
{
	public class ScalingMatrixRowMajor : IMappingMatrix
	{
		public ScalingMatrixRowMajor(int numRows, int numColumns, int[] nonZeroColOfRow, double[] nonZeroValOfRow)
		{
			this.NumRows = numRows;
			this.NumColumns = numColumns;
			this.NonZeroColOfRow = nonZeroColOfRow;
			this.NonZeroValOfRow = nonZeroValOfRow;
		}

		public int[] NonZeroColOfRow { get; }

		public double[] NonZeroValOfRow { get; }

		public int NumColumns { get; }

		public int NumRows { get; }

		public Matrix CopyToFullMatrix()
		{
			var full = Matrix.CreateZero(NumRows, NumColumns);
			for (int i = 0; i < NumRows; ++i)
			{
				full[i, NonZeroColOfRow[i]] = NonZeroValOfRow[i];
			}
			return full;
		}

		public Vector Multiply(Vector vector, bool transpose)
		{
			if (transpose)
			{
				Preconditions.CheckMultiplicationDimensions(NumRows, vector.Length);
				var result = new double[NumColumns];
				for (int i = 0; i < NumRows; ++i)
				{
					result[NonZeroColOfRow[i]] = vector[i] * NonZeroValOfRow[i];
				}
				return Vector.CreateFromArray(result);
			}
			else
			{
				Preconditions.CheckMultiplicationDimensions(NumColumns, vector.Length);
				var result = new double[NumRows];
				for (int i = 0; i < NumRows; ++i)
				{
					result[i] = vector[NonZeroColOfRow[i]] * NonZeroValOfRow[i];
				}
				return Vector.CreateFromArray(result);
			}
		}
	}
}
