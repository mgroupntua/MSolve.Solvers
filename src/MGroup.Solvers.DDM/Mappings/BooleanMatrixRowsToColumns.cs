using System;
using System.Collections.Generic;
using System.Text;
using MGroup.LinearAlgebra.Commons;
using MGroup.LinearAlgebra.Matrices;
using MGroup.LinearAlgebra.Vectors;

namespace MGroup.Solvers.DDM.Mappings
{
	public class BooleanMatrixRowsToColumns : IMappingMatrix
	{
		public BooleanMatrixRowsToColumns(int numRows, int numColumns, int[] rowsToColumns)
		{
			this.NumRows = numRows;
			this.NumColumns = numColumns;
			this.RowsToColumns = rowsToColumns;
		}

		public int[] RowsToColumns { get; }

		public int NumColumns { get; }

		public int NumRows { get; }

		public Matrix CopyToFullMatrix()
		{
			var full = Matrix.CreateZero(NumRows, NumColumns);
			for (int i = 0; i < NumRows; ++i)
			{
				full[i, RowsToColumns[i]] = 1.0;
			}
			return full;
		}

		public Vector Multiply(Vector vector, bool transpose)
		{
			if (transpose)
			{
				Preconditions.CheckMultiplicationDimensions(NumRows, vector.Length);
				var result = new double[NumColumns];
				for (int i = 0; i < NumRows; i++)
				{
					result[RowsToColumns[i]] = vector[i];
				}
				return Vector.CreateFromArray(result);
			}
			else
			{
				Preconditions.CheckMultiplicationDimensions(NumColumns, vector.Length);
				var result = new double[NumRows];
				for (int i = 0; i < NumRows; i++)
				{
					result[i] = vector[RowsToColumns[i]];
				}
				return Vector.CreateFromArray(result);
			}
		}
	}
}
