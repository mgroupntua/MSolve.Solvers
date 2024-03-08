using System;
using System.Collections.Generic;
using System.Text;
using MGroup.LinearAlgebra.Distributed.Overlapping;

namespace MGroup.Solvers.DDM.PSM.Reanalysis
{
	/// <summary>
	/// Will always start the new solution from 0.
	/// </summary>
	public class ZeroInitialSolutionGuess : IInitialSolutionGuessStrategy
	{
		public (DistributedOverlappingVector guess, bool isZero) GuessFirstSolution(
			DistributedOverlappingIndexer currentBoundaryDofIndexer)
		{
			return (new DistributedOverlappingVector(currentBoundaryDofIndexer), true);
		}

		public (DistributedOverlappingVector guess, bool isZero) GuessNextSolution(
			DistributedOverlappingIndexer currentBoundaryDofIndexer, DistributedOverlappingVector previousSolution)
		{
			return (new DistributedOverlappingVector(currentBoundaryDofIndexer), true);
		}
	}
}
