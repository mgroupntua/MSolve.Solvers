using System;
using System.Collections.Generic;
using System.Runtime.CompilerServices;
using System.Text;

namespace MGroup.Solvers.DDM.LinearAlgebraExtensions
{
	public class IntSymMatrixColMajor
	{
		private readonly int order;

		/// <summary>
		/// Packed storage, column major order, upper triangle: 
		/// A[i,j] = data[i + j*(j+1)/2] for 0 &lt;= i &lt;= j &lt; n.
		/// </summary>
		private readonly int[] values;

		private IntSymMatrixColMajor(int order, int[] values)
		{
			this.order = order;
			this.values = values;
			this.NumRows = order;
			this.NumColumns = order;
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
		/// Packed storage, column major order, upper triangle: 
		/// A[i,j] = data[i + j*(j+1)/2] for 0 &lt;= i &lt;= j &lt; n.
		/// </summary>
		public int[] RawValues => values;

		public int this[int rowIdx, int colIdx]
		{
			get
			{
				if (rowIdx <= colIdx)
				{
					return values[Find1DIndex(rowIdx, colIdx)];
				}
				else
				{
					return values[Find1DIndex(colIdx, rowIdx)];
				}
			}

			set
			{
				if (rowIdx <= colIdx)
				{
					values[Find1DIndex(rowIdx, colIdx)] = value;
				}
				else
				{
					values[Find1DIndex(colIdx, rowIdx)] = value;
				}
			}
		}

		/// <summary>
		/// Create a new <see cref="IntSymMatrixColMajor"/> with the specified order and all entries equal to 0.
		/// </summary> 
		/// <param name="order">The number of rows or columns of the matrix.</param>
		public static IntSymMatrixColMajor CreateZero(int order)
		{
			var values = new int[((order + 1) * order) / 2];
			return new IntSymMatrixColMajor(order, values);
		}

		[MethodImpl(MethodImplOptions.AggressiveInlining)]
		internal int Find1DIndex(int i, int j) => i + (j * (j + 1)) / 2;
	}
}
