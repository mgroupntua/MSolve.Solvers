using System;
using System.Collections.Generic;
using System.Text;

namespace MGroup.Solvers.DDM.LinearAlgebraExtensions
{
	public class IntDokRowMajor
	{
		private readonly SortedDictionary<int, int>[] rows;

		public IntDokRowMajor(int numRows, int numColumns, SortedDictionary<int, int>[] columns)
		{
			this.rows = columns;
			this.NumRows = numRows;
			this.NumColumns = numColumns;
		}

		public int NumColumns { get; }

		public int NumRows { get; }

		public static IntDokRowMajor CreateZero(int numRows, int numColumns)
		{
			var rows = new SortedDictionary<int, int>[numRows];
			for (int i = 0; i < numRows; ++i)
			{
				rows[i] = new SortedDictionary<int, int>(); //Initial capacity may be optimized.
			}
			return new IntDokRowMajor(numRows, numColumns, rows);
		}

		public int this[int rowIdx, int colIdx]
		{
			get
			{
				if (rows[rowIdx].TryGetValue(colIdx, out int val))
				{
					return val;
				}
				else return 0;
			}
			set //not thread safe
			{
				rows[rowIdx][colIdx] = value;
			}
		}

		public (int[] values, int[] colIndices, int[] rowOffsets) BuildCsrArrays()
		{
			int[] rowOffsets = new int[NumRows + 1];
			int nnz = 0;
			for (int i = 0; i < NumRows; ++i)
			{
				rowOffsets[i] = nnz;
				nnz += rows[i].Count;
			}
			if (nnz == 0) throw new Exception("Cannot build symmetric CSC arrays from a DOK with nnz = 0.");
			rowOffsets[NumRows] = nnz; //The last CSR entry has colOffset = nnz.

			int[] colIndices = new int[nnz];
			int[] values = new int[nnz];
			int counter = 0;
			for (int i = 0; i < NumRows; ++i)
			{
				// These are already ordered by the SortedDictionary
				foreach (var colValPair in rows[i])
				{
					colIndices[counter] = colValPair.Key;
					values[counter] = colValPair.Value;
					++counter;
				}
			}
			return (values, colIndices, rowOffsets);
		}
	}
}
