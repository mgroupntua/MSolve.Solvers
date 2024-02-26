using System;
using System.Collections.Generic;
using System.Text;

namespace MGroup.Solvers.DDM.LinearAlgebraExtensions
{
	public class IntDokSymColMajor
	{
		private readonly int order;
		private readonly SortedDictionary<int, int>[] columns;

		public IntDokSymColMajor(int order, SortedDictionary<int, int>[] columns)
		{
			this.order = order;
			this.columns = columns;

			this.NumRows = order;
			this.NumColumns = order;
		}

		public int NumColumns { get; }

		public int NumRows { get; }

		public static IntDokSymColMajor CreateZero(int order)
		{
			var columns = new SortedDictionary<int, int>[order];
			for (int j = 0; j < order; ++j)
			{
				columns[j] = new SortedDictionary<int, int>(); //Initial capacity may be optimized.
			}
			return new IntDokSymColMajor(order, columns);
		}

		public int this[int rowIdx, int colIdx]
		{
			get
			{
				if (rowIdx > colIdx)
				{
					int swap = rowIdx;
					rowIdx = colIdx;
					colIdx = swap;
				}

				if (columns[colIdx].TryGetValue(rowIdx, out int val))
				{
					return val;
				}
				else
				{
					return 0;
				}
			}

			set //not thread safe
			{
				if (rowIdx > colIdx)
				{
					int swap = rowIdx;
					rowIdx = colIdx;
					colIdx = swap;
				}

				columns[colIdx][rowIdx] = value;
			}
		}

		public (int[] values, int[] rowIndices, int[] colOffsets) BuildCscArrays()
		{
			int[] colOffsets = new int[order + 1];
			int nnz = 0;
			for (int j = 0; j < order; ++j)
			{
				colOffsets[j] = nnz;
				nnz += columns[j].Count;
			}
			if (nnz == 0) throw new Exception("Cannot build symmetric CSC arrays from a DOK with nnz = 0.");
			colOffsets[order] = nnz; //The last CSC entry has colOffset = nnz.

			int[] rowIndices = new int[nnz];
			int[] values = new int[nnz];
			int counter = 0;
			for (int j = 0; j < order; ++j)
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
