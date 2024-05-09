namespace MGroup.Solvers.DDM.PSM.InterfaceProblem
{
	using MGroup.LinearAlgebra.Distributed.Overlapping;

	public interface IPsmInterfaceProblemMatrix
	{
		DistributedOverlappingTransformation Matrix { get; }

		void Calculate(DistributedOverlappingIndexer indexer);

		double[] ExtractDiagonal(int subdomainID);
	}
}
