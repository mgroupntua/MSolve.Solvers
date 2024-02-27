using System;
using System.Collections.Generic;
using System.Text;
using MGroup.LinearAlgebra.Commons;
using MGroup.LinearAlgebra.Matrices;
using MGroup.LinearAlgebra.Vectors;

namespace MGroup.Solvers.DDM.Mappings
{
	public class MappingMatrixN : IMappingMatrix
	{
		public MappingMatrixN(int numRows, int numColumns, Dictionary<int, int> rowsToColumns)
		{
			this.RowsToColumns = rowsToColumns;
			this.NumRows = numRows;
			this.NumColumns = numColumns;
		}

		public Dictionary<int, int> RowsToColumns { get; }

		public int NumColumns { get; }

		public int NumRows { get; }

		public Matrix CopyToFullMatrix()
		{
			var full = Matrix.CreateZero(NumRows, NumColumns);
			foreach (var rowCol in RowsToColumns)
			{
				full[rowCol.Key, rowCol.Value] = 1.0;
			}
			return full;
		}

		public Vector Multiply(Vector vector, bool transpose)
		{
			if (transpose)
			{
				Preconditions.CheckMultiplicationDimensions(NumRows, vector.Length);
				var result = new double[NumColumns];
				foreach (var rowCol in RowsToColumns)
				{
					result[rowCol.Value] = vector[rowCol.Key];
				}
				return Vector.CreateFromArray(result);
			}
			else
			{
				Preconditions.CheckMultiplicationDimensions(NumColumns, vector.Length);
				var result = new double[NumRows];
				foreach (var rowCol in RowsToColumns)
				{
					result[rowCol.Key] = vector[rowCol.Value];
				}
				return Vector.CreateFromArray(result);
			}
		}
	}
}
