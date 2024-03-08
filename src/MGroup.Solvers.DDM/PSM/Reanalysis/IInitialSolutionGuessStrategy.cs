using System;
using System.Collections.Generic;
using System.Text;
using MGroup.LinearAlgebra.Distributed.Overlapping;
using MGroup.Solvers.DDM.LinearSystem;

namespace MGroup.Solvers.DDM.PSM.Reanalysis
{
	public interface IInitialSolutionGuessStrategy
	{
		(DistributedOverlappingVector guess, bool isZero) GuessFirstSolution(
			DistributedOverlappingIndexer currentBoundaryDofIndexer);

		(DistributedOverlappingVector guess, bool isZero) GuessNextSolution(
			DistributedOverlappingIndexer currentBoundaryDofIndexer, DistributedOverlappingVector previousSolution);
	}
}
