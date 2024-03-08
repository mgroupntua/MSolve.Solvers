using MGroup.LinearAlgebra.Distributed.Overlapping;

namespace MGroup.Solvers.DDM.FetiDP.InterfaceProblem
{
	public interface IFetiDPInterfaceProblemVectors
	{
		DistributedOverlappingVector InterfaceProblemRhs { get; }

		DistributedOverlappingVector InterfaceProblemSolution { get; set; }

		void CalcInterfaceRhsVector(DistributedOverlappingIndexer lagrangeVectorIndexer);

		void Clear();
	}
}
