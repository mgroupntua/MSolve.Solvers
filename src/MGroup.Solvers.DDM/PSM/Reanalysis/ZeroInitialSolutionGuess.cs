namespace MGroup.Solvers.DDM.PSM.Reanalysis
{
	using MGroup.LinearAlgebra.Distributed.Overlapping;

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
