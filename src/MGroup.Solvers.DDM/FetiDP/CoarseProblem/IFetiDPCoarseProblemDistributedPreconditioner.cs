using MGroup.LinearAlgebra.Distributed.IterativeMethods.Preconditioning;
using MGroup.LinearAlgebra.Distributed.Overlapping;

namespace MGroup.Solvers.DDM.FetiDP.CoarseProblem
{
	public interface IFetiDPCoarseProblemDistributedPreconditioner
	{
		IPreconditioner Preconditioner { get; }

		void Calculate(DistributedOverlappingIndexer cornerDofIndexer);
	}
}
