namespace MGroup.Solvers.DDM.FetiDP.InterfaceProblem
{
	using MGroup.LinearAlgebra.Distributed.Overlapping;

	public interface IFetiDPInterfaceProblemVectors
	{
		DistributedOverlappingVector InterfaceProblemRhs { get; }

		DistributedOverlappingVector InterfaceProblemSolution { get; set; }

		void CalcInterfaceRhsVector(DistributedOverlappingIndexer lagrangeVectorIndexer);

		void Clear();
	}
}
