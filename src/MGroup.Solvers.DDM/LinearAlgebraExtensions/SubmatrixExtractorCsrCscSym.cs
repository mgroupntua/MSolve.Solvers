using System;
using System.Collections.Generic;
using System.Text;
using MGroup.LinearAlgebra.Matrices;

namespace MGroup.Solvers.DDM.LinearAlgebraExtensions
{
	public class SubmatrixExtractorCsrCscSym : SubmatrixExtractorCscSymBase
	{
		public SubmatrixExtractorCsrCscSym()
		{
		}

		public CsrMatrix Submatrix00 { get; private set; }

		public override void Clear()
		{
			base.Clear();
			Submatrix00 = null;
		}


		public void ExtractSubmatrices(SymmetricCscMatrix originalMatrix, int[] indicesGroup0, int[] indicesGroup1)
		{
			if (this.originalMatrix == null)
			{
				ExtractMaps(originalMatrix, indicesGroup0, indicesGroup1);
				this.originalMatrix = originalMatrix;
			}
			else
			{
				if (this.originalMatrix != originalMatrix)
				{
					throw new InvalidOperationException(
						"The submatrix pattern stored corresponds to a different original matrix. Clear this object first");
				}
			}

			CopyValuesArray(originalMatrix.RawValues, Submatrix00.RawValues, map00);
			CopyValuesArray(originalMatrix.RawValues, Submatrix01.RawValues, map01);
			CopyValuesArray(originalMatrix.RawValues, Submatrix11.RawValues, map11);
		}

		private void ExtractMaps(SymmetricCscMatrix originalMatrix, int[] indicesGroup0, int[] indicesGroup1)
		{
			int n0 = indicesGroup0.Length;
			int n1 = indicesGroup1.Length;
			var submatrix00 = IntDokRowMajor.CreateZero(n0, n0);
			var submatrix01 = IntDokRowMajor.CreateZero(n0, n1);
			var submatrix11 = IntDokSymColMajor.CreateZero(n1);

			// Original matrix indices to submatrix indices. Group 0 indices i0 are stored as: -i0-1. 
			// Group 1 indices are stored as they are
			int[] originalToSubIndices = MapOriginalToSubmatrixIndices(originalMatrix.NumRows, indicesGroup0, indicesGroup1);

			// Iterate the non zero values array of the original sparse matrix
			for (int j = 0; j < originalMatrix.NumColumns; ++j)
			{
				int start = originalMatrix.RawColOffsets[j];
				int end = originalMatrix.RawColOffsets[j + 1];
				int subJ = originalToSubIndices[j];
				if (subJ < 0) // j belongs to group 0 indices
				{
					subJ = -subJ - 1; // Mapping to [0, oo) for group 0 indices
					for (int t = start; t < end; ++t)
					{
						int i = originalMatrix.RawRowIndices[t];
						int subI = originalToSubIndices[i];
						if (subI < 0) // (i, j) belongs to submatrix A00
						{
							subI = -subI- 1; // Mapping to [0, oo) for group 0 indices
							submatrix00[subI, subJ] = t;
							submatrix00[subJ, subI] = t;
						}
						else // (i, j) belongs to submatrix A10 or equivalently (j, i) belongs to submatrix A01
						{
							submatrix01[subJ, subI] = t;
						}
					}
				}
				else // j belongs to group 1 indices
				{
					for (int t = start; t < end; ++t)
					{
						int i = originalMatrix.RawRowIndices[t];
						int subI = originalToSubIndices[i];
						if (subI < 0) // (i, j) belongs to submatrix A01
						{
							subI = -subI- 1; // Mapping to [0, oo) for group 0 indices
							submatrix01[subI, subJ] = t;
						}
						else // (i, j) belongs to submatrix A11
						{
							submatrix11[subI, subJ] = t;
						}
					}
				}
			}

			// Finalize the data structures required to represent the submatrices
			// A00 CSR
			int[] colIndices00, rowOffsets00;
			(this.map00, colIndices00, rowOffsets00) = submatrix00.BuildCsrArrays();
			this.Submatrix00 = CsrMatrix.CreateFromArrays(
				n0, n0, new double[this.map00.Length], colIndices00, rowOffsets00, false);

			// A01 CSR
			int[] colIndices01, rowOffsets01;
			(this.map01, colIndices01, rowOffsets01) = submatrix01.BuildCsrArrays();
			this.Submatrix01 = CsrMatrix.CreateFromArrays(
				n0, n1, new double[this.map01.Length], colIndices01, rowOffsets01, false);

			// A11 CSC upper triangle
			int[] rowIndices11, colOffsets11;
			(this.map11, rowIndices11, colOffsets11) = submatrix11.BuildCscArrays();
			this.Submatrix11 = SymmetricCscMatrix.CreateFromArrays(
				n1, new double[this.map11.Length], rowIndices11, colOffsets11, false);
		}
	}
}
