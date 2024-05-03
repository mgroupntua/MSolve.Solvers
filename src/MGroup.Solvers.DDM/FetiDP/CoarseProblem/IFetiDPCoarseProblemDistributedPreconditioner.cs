namespace MGroup.Solvers.DDM.FetiDP.CoarseProblem
{
	using MGroup.LinearAlgebra.Distributed.IterativeMethods.Preconditioning;
	using MGroup.LinearAlgebra.Distributed.Overlapping;

	public interface IFetiDPCoarseProblemDistributedPreconditioner
	{
		IPreconditioner Preconditioner { get; }

		void Calculate(DistributedOverlappingIndexer cornerDofIndexer);
	}
}
