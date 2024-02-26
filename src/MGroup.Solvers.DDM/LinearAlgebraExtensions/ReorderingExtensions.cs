using System;
using System.Collections.Generic;
using System.Runtime.CompilerServices;
using System.Text;
using CSparse;
using CSparse.Double;
using CSparse.Ordering;
using MGroup.LinearAlgebra.Reordering;

namespace MGroup.Solvers.DDM.LinearAlgebraExtensions
{
	public static class ReorderingExtensions
	{
		/// <summary>
		/// See <see cref="IReorderingAlgorithm.FindPermutation(SparsityPatternSymmetric)"/>
		/// </summary>
		/// <remarks>The returned permutation is new-to-old.</remarks>
		public static (int[] permutation, bool oldToNew) FindPermutation(this OrderingAmdCSparseNet reordering, int order,
			int[] cscRowIndices, int[] cscColOffsets)
		{
			var dummyCscValues = new double[cscRowIndices.Length]; //TODO: too expensive 
			var matrixCSparse = new SparseMatrix(order, order, dummyCscValues, cscRowIndices, cscColOffsets);
			int[] permutation = AMD.Generate<double>(matrixCSparse, ColumnOrdering.MinimumDegreeAtPlusA);

			// It is possible that CSparse.NET AMD algorithm returns more entries than the matrix order (so far I have found 1 
			// extra). In that case, make sure the first ones are valid and return only them.
			if (permutation.Length > order)
			{
				for (int i = order; i < permutation.Length; ++i)
				{
					if (permutation[i] < order) throw new Exception(
						"Something went wrong during AMD. The permutation vector has more entries than the matrix order.");
				}
				var permutationCorrected = new int[order];
				Array.Copy(permutation, permutationCorrected, order);
				return (permutationCorrected, false);
			}
			else return (permutation, false);
		}
	}
}
