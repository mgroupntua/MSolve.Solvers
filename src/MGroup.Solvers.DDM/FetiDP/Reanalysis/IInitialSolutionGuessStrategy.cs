using MGroup.LinearAlgebra.Distributed.Overlapping;

namespace MGroup.Solvers.DDM.FetiDP.Reanalysis
{
	//TODO: I could reuse the one of PSM
	public interface IInitialSolutionGuessStrategy
	{
		(DistributedOverlappingVector guess, bool isZero) GuessFirstSolution(
			DistributedOverlappingIndexer currentLagrangeVectorIndexer);

		(DistributedOverlappingVector guess, bool isZero) GuessNextSolution(
			DistributedOverlappingIndexer currentLagrangeVectorIndexer, DistributedOverlappingVector previousSolution);
	}
}
