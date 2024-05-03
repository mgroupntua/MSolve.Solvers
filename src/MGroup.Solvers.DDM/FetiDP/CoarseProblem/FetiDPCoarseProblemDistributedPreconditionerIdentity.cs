namespace MGroup.Solvers.DDM.FetiDP.CoarseProblem
{
	using MGroup.LinearAlgebra.Distributed.IterativeMethods.Preconditioning;
	using MGroup.LinearAlgebra.Distributed.Overlapping;

	public class FetiDPCoarseProblemDistributedPreconditionerIdentity : IFetiDPCoarseProblemDistributedPreconditioner
	{
		public IPreconditioner Preconditioner { get; } = new IdentityPreconditioner();

		public void Calculate(DistributedOverlappingIndexer cornerDofIndexer) { }
	}
}
