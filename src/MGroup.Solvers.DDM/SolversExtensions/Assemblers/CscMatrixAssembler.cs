//TODO: Merge this with the assembler used for element -> subdomain map-reductions.
namespace MGroup.Solvers.DDM.SolversExtensions.Assemblers
{
	using System.Collections.Generic;
	using System.Diagnostics;

	using MGroup.LinearAlgebra.Matrices;
	using MGroup.LinearAlgebra.Matrices.Builders;

	public class CscMatrixAssembler
	{
		private const string name = "CscMatrixAssembler"; // for error messages
		private readonly bool isSymmetric;
		private readonly bool sortColsOfEachRow;

		bool isIndexerCached = false;
		private int[] cachedRowIndices, cachedColOffsets;

		/// <summary>
		/// 
		/// </summary>
		/// <param name="sortColsOfEachRow">
		/// Sorting the columns of each row in the CSC storage format may increase performance of the factorization and 
		/// back/forward substitutions. It is recommended to set it to true.
		/// </param>
		public CscMatrixAssembler(bool isSymmetric, bool sortColsOfEachRow = true)
		{
			this.isSymmetric = isSymmetric;
			this.sortColsOfEachRow = sortColsOfEachRow;
		}

		public CscMatrix BuildGlobalMatrix(int numGlobalDofs,
			IDictionary<int, int[]> localToGlobalMaps, IDictionary<int, IMatrix> localMatrices)
		{
			var globalMatrix = DokColMajor.CreateEmpty(numGlobalDofs, numGlobalDofs);
			foreach (var s in localToGlobalMaps.Keys)
			{
				var globalDofs = localToGlobalMaps[s];
				var localDofs = Utilities.Range(0, globalDofs.Length); //TODO: Create a DokSymmetric.AddSubmatrixSymmetric() overload that accepts a single mapping array
				if (isSymmetric)
				{
					globalMatrix.AddSubmatrixSymmetric(localMatrices[s], localDofs, globalDofs);
				}
				else
				{
					globalMatrix.AddSubmatrix(localMatrices[s], localDofs, globalDofs, localDofs, globalDofs);
				}
				globalMatrix.AddSubmatrixSymmetric(localMatrices[s], localDofs, globalDofs);
			}

			(var values, var rowIndices, var colOffsets) = globalMatrix.BuildCscArrays(sortColsOfEachRow);
			if (!isIndexerCached)
			{
				cachedRowIndices = rowIndices;
				cachedColOffsets = colOffsets;
				isIndexerCached = true;
			}
			else
			{
				Debug.Assert(Utilities.AreEqual(cachedRowIndices, rowIndices));
				Debug.Assert(Utilities.AreEqual(cachedColOffsets, colOffsets));
			}

			return CscMatrix.CreateFromArrays(numGlobalDofs, numGlobalDofs, values, cachedRowIndices, cachedColOffsets, false);
		}

		public void HandleDofOrderingWasModified()
		{
			//TODO: perhaps the indexer should be disposed altogether. Then again it could be in use by other matrices.
			cachedRowIndices = null;
			cachedColOffsets = null;
			isIndexerCached = false;
		}
	}
}
