using System.Collections.Generic;
using System.Diagnostics;
using MGroup.LinearAlgebra.Matrices;
using MGroup.LinearAlgebra.Matrices.Builders;

//TODO: Merge this with the assembler used for element -> subdomain map-reductions.
namespace MGroup.Solvers.DDM.AssemblerExtensions
{
	public class SymmetricCscMatrixAssembler
	{
		private const string name = "SymmetricCscMatrixAssembler"; // for error messages
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
		public SymmetricCscMatrixAssembler(bool sortColsOfEachRow = true)
		{
			this.sortColsOfEachRow = sortColsOfEachRow;
		}

		public SymmetricCscMatrix BuildGlobalMatrix(int numGlobalDofs,
			IDictionary<int, int[]> localToGlobalMaps, IDictionary<int, IMatrix> localMatrices)
		{
			var globalMatrix = DokSymmetric.CreateEmpty(numGlobalDofs);
			foreach (int s in localToGlobalMaps.Keys)
			{
				int[] globalDofs = localToGlobalMaps[s];
				int[] localDofs = Utilities.Range(0, globalDofs.Length); //TODO: Create a DokSymmetric.AddSubmatrixSymmetric() overload that accepts a single mapping array
				globalMatrix.AddSubmatrixSymmetric(localMatrices[s], localDofs, globalDofs);
			}

			(double[] values, int[] rowIndices, int[] colOffsets) = globalMatrix.BuildSymmetricCscArrays(sortColsOfEachRow);
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


			return SymmetricCscMatrix.CreateFromArrays(numGlobalDofs, values, cachedRowIndices, cachedColOffsets, false);
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
