using MGroup.LinearAlgebra.Distributed.Overlapping;
using MGroup.MSolve.Solution.LinearSystem;

namespace MGroup.Solvers.DDM.FetiDP.InterfaceProblem
{
	public interface IFetiDPInterfaceProblemMatrix : ILinearTransformation
	{
		void Calculate(DistributedOverlappingIndexer lagrangeVectorIndexer);
	}
}
