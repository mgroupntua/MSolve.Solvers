namespace MGroup.Solvers.DDM.PSM.InterfaceProblem
{
	using MGroup.LinearAlgebra.Distributed.Overlapping;

	public interface IPsmInterfaceProblemVectors
	{
		DistributedOverlappingVector InterfaceProblemRhs { get; }

		DistributedOverlappingVector InterfaceProblemSolution { get; set; }

		// globalF = sum {Lb[s]^T * (fb[s] - Kbi[s] * inv(Kii[s]) * fi[s]) }
		void CalcInterfaceRhsVector(DistributedOverlappingIndexer indexer);

		void Clear();

	}
}
