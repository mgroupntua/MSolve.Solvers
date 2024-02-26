using System;
using System.Collections.Generic;
using System.Text;
using MGroup.LinearAlgebra.Matrices;
using MGroup.LinearAlgebra.Triangulation;
using MGroup.LinearAlgebra.Vectors;

namespace MGroup.Solvers.DDM.LinearAlgebraExtensions
{
	public static class SchurComplementPckCsrCscSym
	{
		/// <summary>
		/// Calculates the Schur complement of A/A22 = S = A11 - A21^T * inv(A22) * A21, where M = [A11 A12; A12^T A22].
		/// This method constructs inv(A22) * A12^T one column at a time and uses that column to calculate the superdiagonal
		/// entries of the corresponding column of A12 * inv(A22) * A12^T.
		/// </summary>
		public static SymmetricMatrix CalcSchurComplement(SymmetricMatrix A11, CsrMatrix A12, ITriangulation inverseA22)
		{
			var S = SymmetricMatrix.CreateZero(A11.Order);
			CalcSchurComplement(A11, A12, inverseA22, S);
			return S;
		}

		/// <summary>
		/// Calculates the Schur complement of A/A22 = S = A11 - A21^T * inv(A22) * A21, where M = [A11 A12; A12^T A22].
		/// This method constructs inv(A22) * A12^T one column at a time and uses that column to calculate the superdiagonal
		/// entries of the corresponding column of A12 * inv(A22) * A12^T.
		/// </summary>
		public static void CalcSchurComplement(SymmetricMatrix A11, CsrMatrix A12, ITriangulation inverseA22, 
			SymmetricMatrix result)
		{ //TODO: Unfortunately this cannot take advantage of MKL for CSR * matrix or even CSR * vector.
			// Column j of A21 = row j of A12
			for (int j = 0; j < A12.NumRows; ++j)
			{
				// column j of (inv(A22) * A21) = inv(A22) * column j of A21
				Vector colA21 = A12.GetRow(j);
				double[] colInvA22A21 = inverseA22.SolveLinearSystem(colA21).RawData;

				// column j of (A12 * inv(A22) * A21) = A12 * column j of (inv(A22) * A21)
				// However we only need the superdiagonal part of this column. 
				// Thus we only multiply the rows i of A12 with i <= j. 
				for (int i = 0; i <= j; ++i)
				{
					// Perform the subtraction S = A11 - (A21^T * inv(A22) * A21) for the current (i, j)
					double dot = A12.MultiplyRowTimesVector(i, colInvA22A21);
					int indexS = i + (j * (j + 1)) / 2;
					result.RawData[indexS] = A11.RawData[indexS] - dot;
				}
			}
		}
	}
}
