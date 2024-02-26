using System;
using System.Collections.Generic;
using System.Text;
using MGroup.LinearAlgebra.Matrices;
using MGroup.LinearAlgebra.Matrices.Builders;

namespace MGroup.Solvers.DDM.LinearAlgebraExtensions
{
	public static class DokExtensions
	{
		public static CsrMatrix GetSubmatrixCsr(this DokRowMajor dok, int[] rows, int[] cols)
		{
			var submatrix = DokRowMajor.CreateEmpty(rows.Length, cols.Length);
			for (int i = 0; i < rows.Length; i++)
			{
				int I = rows[i];
				for (int j = 0; j < cols.Length; j++)
				{
					int J = cols[j];
					submatrix[i, j] = dok[I, J];
				}
			}
			return submatrix.BuildCsrMatrix(true);
		}

		public static CscMatrix GetSubmatrixCsc(this DokRowMajor dok, int[] rows, int[] cols)
		{
			var submatrix = DokColMajor.CreateEmpty(rows.Length, cols.Length);
			for (int i = 0; i < rows.Length; i++)
			{
				int I = rows[i];
				for (int j = 0; j < cols.Length; j++)
				{
					int J = cols[j];
					submatrix[i, j] = dok[I, J];
				}
			}
			return submatrix.BuildCscMatrix(true);
		}

		public static CsrMatrix GetSubmatrixCsr(this IIndexable2D originalMatrix, int[] rows, int[] cols)
		{
			var submatrix = DokRowMajor.CreateEmpty(rows.Length, cols.Length);
			for (int i = 0; i < rows.Length; i++)
			{
				int I = rows[i];
				for (int j = 0; j < cols.Length; j++)
				{
					int J = cols[j];
					submatrix[i, j] = originalMatrix[I, J];
				}
			}
			return submatrix.BuildCsrMatrix(true);
		}

		public static CscMatrix GetSubmatrixCsc(this IIndexable2D originalMatrix, int[] rows, int[] cols)
		{
			var submatrix = DokColMajor.CreateEmpty(rows.Length, cols.Length);
			for (int i = 0; i < rows.Length; i++)
			{
				int I = rows[i];
				for (int j = 0; j < cols.Length; j++)
				{
					int J = cols[j];
					submatrix[i, j] = originalMatrix[I, J];
				}
			}
			return submatrix.BuildCscMatrix(true);
		}
	}
}
