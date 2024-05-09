namespace MGroup.Solvers.DDM.PSM.Reanalysis
{
	using MGroup.LinearAlgebra.Distributed.Overlapping;

	public interface IInitialSolutionGuessStrategy
	{
		(DistributedOverlappingVector guess, bool isZero) GuessFirstSolution(
			DistributedOverlappingIndexer currentBoundaryDofIndexer);

		(DistributedOverlappingVector guess, bool isZero) GuessNextSolution(
			DistributedOverlappingIndexer currentBoundaryDofIndexer, DistributedOverlappingVector previousSolution);
	}
}
