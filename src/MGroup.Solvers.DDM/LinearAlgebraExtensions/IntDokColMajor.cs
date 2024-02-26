using System;
using System.Collections.Generic;
using System.Text;

namespace MGroup.Solvers.DDM.LinearAlgebraExtensions
{
	public class IntDokColMajor
	{
		private readonly SortedDictionary<int, int>[] columns;

		public IntDokColMajor(int numRows, int numColumns, SortedDictionary<int, int>[] columns)
		{
			this.columns = columns;
			this.NumRows = numRows;
			this.NumColumns = numColumns;
		}

		public int NumColumns { get; }

		public int NumRows { get; }

		public static IntDokColMajor CreateZero(int numRows, int numColumns)
		{
			var columns = new SortedDictionary<int, int>[numColumns];
			for (int j = 0; j < numColumns; ++j)
			{
				columns[j] = new SortedDictionary<int, int>(); //Initial capacity may be optimized.
			}
			return new IntDokColMajor(numRows, numColumns, columns);
		}

		public int this[int rowIdx, int colIdx]
		{
			get
			{
				if (columns[colIdx].TryGetValue(rowIdx, out int val))
				{
					return val;
				}
				else return 0;
			}
			set //not thread safe
			{
				columns[colIdx][rowIdx] = value;
			}
		}

		public (int[] values, int[] rowIndices, int[] colOffsets) BuildCscArrays()
		{
			int[] colOffsets = new int[NumColumns + 1];
			int nnz = 0;
			for (int j = 0; j < NumColumns; ++j)
			{
				colOffsets[j] = nnz;
				nnz += columns[j].Count;
			}
			if (nnz == 0) throw new Exception("Cannot build symmetric CSC arrays from a DOK with nnz = 0.");
			colOffsets[NumColumns] = nnz; //The last CSC entry has colOffset = nnz.

			int[] rowIndices = new int[nnz];
			int[] values = new int[nnz];
			int counter = 0;
			for (int j = 0; j < NumColumns; ++j)
			{
				// These are already ordered by the SortedDictionary
				foreach (var rowValPair in columns[j])
				{
					rowIndices[counter] = rowValPair.Key;
					values[counter] = rowValPair.Value;
					++counter;
				}
			}
			return (values, rowIndices, colOffsets);
		}
	}
}
