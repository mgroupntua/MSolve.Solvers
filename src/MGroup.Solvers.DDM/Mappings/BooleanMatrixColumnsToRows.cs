using System;
using System.Collections.Generic;
using System.Text;
using MGroup.LinearAlgebra.Commons;
using MGroup.LinearAlgebra.Matrices;
using MGroup.LinearAlgebra.Vectors;

namespace MGroup.Solvers.DDM.Mappings
{
	public class BooleanMatrixColumnsToRows : IMappingMatrix
	{
		public BooleanMatrixColumnsToRows(int numRows, int numColumns, int[] columnsToRows)
		{
			this.NumRows = numRows;
			this.NumColumns = numColumns;
			this.ColumnsToRows = columnsToRows;
		}

		public int[] ColumnsToRows { get; }

		public int NumColumns { get; }

		public int NumRows { get; }

		public Matrix CopyToFullMatrix()
		{
			var full = Matrix.CreateZero(NumRows, NumColumns);
			for (int j = 0; j < NumColumns; ++j)
			{
				full[ColumnsToRows[j], j] = 1.0;
			}
			return full;
		}

		public Vector Multiply(Vector vector, bool transpose)
		{
			if (transpose)
			{
				Preconditions.CheckMultiplicationDimensions(NumRows, vector.Length);
				var result = new double[NumColumns];
				for (int j = 0; j < NumColumns; ++j)
				{
					result[j] = vector[ColumnsToRows[j]];
				}
				return Vector.CreateFromArray(result);
			}
			else
			{
				Preconditions.CheckMultiplicationDimensions(NumColumns, vector.Length);
				var result = new double[NumRows];
				for (int j = 0; j < NumColumns; ++j)
				{
					result[ColumnsToRows[j]] = vector[j];
				}
				return Vector.CreateFromArray(result);
			}
		}
	}
}
