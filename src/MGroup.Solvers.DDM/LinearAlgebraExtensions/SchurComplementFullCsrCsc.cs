using System;
using System.Collections.Generic;
using System.Text;
using MGroup.LinearAlgebra.Matrices;
using MGroup.LinearAlgebra.Triangulation;
using MGroup.LinearAlgebra.Vectors;

namespace MGroup.Solvers.DDM.LinearAlgebraExtensions
{
	public static class SchurComplementFullCsrCsc
	{
		/// <summary>
		/// Calculates the Schur complement of A/A22 = S = A11 - A12 * inv(A22) * A21, where M = [A11 A12; A21 A22].
		/// This method constructs inv(A22) * A21 one column at a time and uses that column to calculate the superdiagonal
		/// entries of the corresponding column of A12 * inv(A22) * A21.
		/// </summary>
		public static Matrix CalcSchurComplement(Matrix A11, CsrMatrix A12, CsrMatrix A21, ITriangulation inverseA22)
		{
			var S = Matrix.CreateZero(A11.NumRows, A11.NumColumns);
			CalcSchurComplement(A11, A12, A21, inverseA22, S);
			return S;
		}

		///// <summary>
		///// Calculates the Schur complement of A/A22 = S = A11 - A12 * inv(A22) * A21, where M = [A11 A12; A21 A22].
		///// This method constructs inv(A22) * A21 one column at a time and uses that column to calculate the superdiagonal
		///// entries of the corresponding column of A12 * inv(A22) * A21.
		///// The rows of each CSR matrix must contain columns in ascending order.
		///// </summary>
		//public static void CalcSchurComplement(Matrix A11, CsrMatrix A12, CsrMatrix A21, ITriangulation inverseA22, 
		//	Matrix result)
		//{ 
		//  //TODO: Assert that rows of CSR matrices are in ascending order

		//	result.CopyFrom(A11); // TODO: Needs option for directly overwritting A11

		//	// Use an auxilliary array to reduce misses when searching column j of each row
		//	var rowStarts = new int[A21.RawRowOffsets.Length];
		//	Array.Copy(A21.RawRowOffsets, rowStarts, rowStarts.Length);

		//	// Column j of A21
		//	for (int j = 0; j < A21.NumColumns; ++j)
		//	{
		//		// column j of (inv(A22) * A21) = inv(A22) * column j of A21
		//		Vector colA21 = A21.GetColumn(j);
		//		Vector colInvA22A21 = inverseA22.SolveLinearSystem(colA21);

		//		// column j of (A12 * inv(A22) * A21) = A12 * column j of (inv(A22) * A21)
		//		Vector colS = A12.Multiply(colInvA22A21);
		//		result.SubtractColumn(j, colS);
		//	}
		//}

		/// <summary>
		/// Calculates the Schur complement of A/A22 = S = A11 - A12 * inv(A22) * A21, where M = [A11 A12; A21 A22].
		/// This method constructs inv(A22) * A21 one column at a time and uses that column to calculate the superdiagonal
		/// entries of the corresponding column of A12 * inv(A22) * A21.
		/// The rows of each CSR matrix must contain columns in ascending order.
		/// </summary>
		public static void CalcSchurComplement(Matrix A11, CsrMatrix A12, CsrMatrix A21, ITriangulation inverseA22,
			Matrix result)
		{
			//TODO: Assert that rows of CSR matrices are in ascending order

			// Use an auxilliary arrays to reduce misses when searching column j of each row
			var firstCols = new int[A21.NumRows]; // Start from col 0 for each row
			var rowStarts = new int[A21.NumRows];
			Array.Copy(A21.RawRowOffsets, rowStarts, rowStarts.Length);

			// Column j of A21
			for (int j = 0; j < A21.NumColumns; ++j)
			{
				// column j of (inv(A22) * A21) = inv(A22) * column j of A21
				Vector colA21 = GetCsrColumn(A21, j, rowStarts, firstCols);
				Vector colInvA22A21 = inverseA22.SolveLinearSystem(colA21);

				// column j of (A12 * inv(A22) * A21) = A12 * column j of (inv(A22) * A21)
				Vector colS = A12.Multiply(colInvA22A21);
				result.SetSubcolumn(j, colS);
			}

			// S = A11 - S
			result.LinearCombinationIntoThis(-1, A11, +1);
		}

		/// <summary>
		/// Efficient implementation that 1) assumes CSR.colIndices is sorted per row, 2) uses and modifies an auxilliary array 
		/// <paramref name="rowStarts"/> to hold the starting index for each row to search for <paramref name="colIdx"/>
		/// </summary>
		/// <param name="csr"></param>
		/// <param name="colIdx"></param>
		/// <param name="rowStarts">Will be modified.</param>
		/// <param name="firstCols">Will be modified.</param>
		private static Vector GetCsrColumn(CsrMatrix csr, int colIdx, int[] rowStarts, int[] firstCols)
		{
			double[] values = csr.RawValues;
			int[] colIndices = csr.RawColIndices;
			int[] rowOffsets = csr.RawRowOffsets;

			var result = new double[csr.NumRows];
			for (int i = 0; i < csr.NumRows; ++i)
			{
				if (colIdx <= firstCols[i])
				{
					// Start searching from rowStarts[i], which is further along rowOffsets[i]. The missing entries were 
					// already searched by previous columns, meaning that they are less than the new col index. 
					int start = rowStarts[i]; //inclusive
					int end = rowOffsets[i + 1]; //exclusive
					for (int k = start; k < end; ++k)
					{
						if (colIndices[k] == colIdx)
						{
							result[i] = values[k];
							rowStarts[i] = k + 1; // In the next search ignore all entries up to and including this one.
							firstCols[i] = colIdx + 1;
							break;
						}
						else if (colIndices[k] > colIdx) // no need to search pass this col index
						{
							rowStarts[i] = k; // In the next search ignore all entries up to but excluding this one.
											  //TODO: E.g say cols 3 and 7 are nonzero corresponding to k = 31, 32. After searching for col 4, 
											  //		rowStarts[i]=32. Searching for cols 5, 6 will repeat the loop and check and set rowStarts[i]=32.
											  //		It would be nice if that could be avoided.
							firstCols[i] = colIndices[k];
							break;
						}
					}
				}
				// otherwise there is no point in searching colIndices
			}

			return Vector.CreateFromArray(result);
		}
	}
}
