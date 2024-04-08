using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Text;
using MGroup.LinearAlgebra.Matrices;

namespace MGroup.Solvers.DDM.LinearAlgebraExtensions
{
	public static class CsrMatrixExtensions
	{
		//TODO: Implement a SchurComplement with A11 being CSR and replace this method. 
		/// <summary>
		/// Extracts the upper triangle of <paramref name="csr"/>.
		/// </summary>
		/// <param name="csr"></param>
		public static TriangularUpper CopyToUpperPacked(this CsrMatrix csr) 
		{
			Debug.Assert(csr.NumRows == csr.NumColumns);
			var pck = TriangularUpper.CreateZero(csr.NumColumns);
			for (int i = 0; i < csr.NumRows; ++i)
			{
				int rowStart = csr.RawRowOffsets[i];
				int rowEnd = csr.RawRowOffsets[i + 1]; //exclusive
				for (int k = rowStart; k < rowEnd; ++k)
				{
					int j = csr.RawColIndices[k];
					if (j >= i) //TODO: optimizations are possible if 
					{
						pck[i, j] = csr.RawValues[k];
					}
				}
			}
			return pck;
		}
	}
}
