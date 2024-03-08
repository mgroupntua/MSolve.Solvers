using MGroup.LinearAlgebra.Distributed.IterativeMethods.Preconditioning;
using MGroup.LinearAlgebra.Distributed.Overlapping;

namespace MGroup.Solvers.DDM.FetiDP.CoarseProblem
{
	public class FetiDPCoarseProblemDistributedPreconditionerIdentity : IFetiDPCoarseProblemDistributedPreconditioner
	{
		public IPreconditioner Preconditioner { get; } = new IdentityPreconditioner();

		public void Calculate(DistributedOverlappingIndexer cornerDofIndexer) { }
	}
}
