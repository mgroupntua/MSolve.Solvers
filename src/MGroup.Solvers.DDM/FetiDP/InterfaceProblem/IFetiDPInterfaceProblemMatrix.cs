namespace MGroup.Solvers.DDM.FetiDP.InterfaceProblem
{
	using MGroup.LinearAlgebra.Distributed.Overlapping;
	using MGroup.MSolve.Solution.LinearSystem;

	public interface IFetiDPInterfaceProblemMatrix : ILinearTransformation
	{
		void Calculate(DistributedOverlappingIndexer lagrangeVectorIndexer);
	}
}
