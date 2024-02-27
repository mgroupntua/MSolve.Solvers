using System.Collections.Generic;
using System.Diagnostics;
using MGroup.LinearAlgebra.Matrices;

//TODO: Merge this with the assembler used for element -> subdomain map-reductions.
namespace MGroup.Solvers.DDM.AssemblerExtensions
{
	public class DenseMatrixAssembler
	{
		public Matrix BuildGlobalMatrix(int numGlobalDofs,
			IDictionary<int, int[]> localToGlobalMaps, IDictionary<int, IMatrix> localMatrices)
		{
			var globalMatrix = Matrix.CreateZero(numGlobalDofs, numGlobalDofs);

			// Process the stiffness of each element
			foreach (int s in localToGlobalMaps.Keys)
			{
				int[] globalDofs = localToGlobalMaps[s];
				int[] localDofs = Utilities.Range(0, globalDofs.Length); //TODO: Create a DokSymmetric.AddSubmatrixSymmetric() overload that accepts a single mapping array
				AddLocalToGlobalMatrix(globalMatrix, localMatrices[s], localDofs, globalDofs);
			}

			return globalMatrix;
		}

		public void HandleDofOrderingWasModified()
		{
		   // Do nothing, since there are no idexing arrays to cache.
		}

		private static void AddLocalToGlobalMatrix(Matrix globalMatrix, IMatrixView localMatrix,
			int[] localIndices, int[] globalIndices)
		{
			Debug.Assert(localMatrix.NumRows == localMatrix.NumColumns);
			Debug.Assert(globalIndices.Length == localIndices.Length);

			int numRelevantRows = localIndices.Length;
			for (int i = 0; i < numRelevantRows; ++i)
			{
				int localRow = localIndices[i];
				int globalRow = globalIndices[i];
				for (int j = 0; j < numRelevantRows; ++j)
				{
					int localCol = localIndices[j];
					int globalCol = globalIndices[j];

					globalMatrix[globalRow, globalCol] += localMatrix[localRow, localCol];
				}
			}
		}
	}
}
