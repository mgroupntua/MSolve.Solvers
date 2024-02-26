using System;
using System.Collections.Generic;
using System.Text;

namespace MGroup.Solvers.DDM.LinearAlgebraExtensions
{
	public class IntMatrixColMajor
	{
		private readonly int[] values;

		private IntMatrixColMajor(int[] values, int numRows, int numColumns)
		{
			this.values = values;
			this.NumRows = numRows;
			this.NumColumns = numColumns;
		}

		/// <summary>
		/// The number of columns of the matrix. 
		/// </summary>
		public int NumColumns { get; }


		/// <summary>
		/// The number of rows of the matrix.
		/// </summary>
		public int NumRows { get; }

		/// <summary>
		/// The internal array that stores the entries of the matrix in column major layout.
		/// It should only be used for passing the raw array to linear algebra libraries.
		/// </summary>
		public int[] RawValues => values;

		public int this[int rowIdx, int colIdx]
		{
			get => values[colIdx * NumRows + rowIdx];
			set => values[colIdx * NumRows + rowIdx] = value;
		}

		/// <summary>
		/// Initializes a new instance of <see cref="IntMatrixColMajor"/> with all entries being equal to 0.
		/// </summary> 
		/// <param name="numRows">The number of rows of the new matrix.</param>
		/// <param name="numColumns">The number of rows of the new matrix.</param>
		public static IntMatrixColMajor CreateZero(int numRows, int numColumns)
		{
			var values = new int[numRows * numColumns];
			return new IntMatrixColMajor(values, numRows, numColumns);
		}
	}
}
