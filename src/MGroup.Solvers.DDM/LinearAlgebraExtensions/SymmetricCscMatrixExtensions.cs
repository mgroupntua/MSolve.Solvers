using System;
using System.Collections.Generic;
using System.Text;
using MGroup.LinearAlgebra.Matrices;
using MGroup.LinearAlgebra.Matrices.Builders;

namespace MGroup.Solvers.DDM.LinearAlgebraExtensions
{
	public static class SymmetricCscMatrixExtensions
	{
		public static CsrMatrix ConvertToCsr(this SymmetricCscMatrix original)
		{
			int order = original.NumColumns;
			double[] values = original.RawValues;
			int[] rowIndices = original.RawRowIndices;
			int[] colOffsets = original.RawColOffsets;

			var dok = DokRowMajor.CreateEmpty(order, order);
			for (int j = 0; j < order; ++j)
			{
				int colStart = colOffsets[j];
				int colEnd = colOffsets[j + 1];
				for (int k = colStart; k < colEnd; ++k)
				{
					int i = rowIndices[k];
					double val = values[k];
					dok[i, j] = val;

					if (i != j)
					{
						dok[j, i] = val;
					}
				}
			}

			return dok.BuildCsrMatrix(true);
		}
	}
}
