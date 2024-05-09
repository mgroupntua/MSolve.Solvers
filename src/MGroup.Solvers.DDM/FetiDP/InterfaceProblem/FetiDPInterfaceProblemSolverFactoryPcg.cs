//TODO: This is basically an adapter because these strategies are not generalized to iterative methods other than PCG.
//		Refactor the iterative methods in LinearAlgebra and remove this kind of adapters.
namespace MGroup.Solvers.DDM.FetiDP.InterfaceProblem
{
	using MGroup.LinearAlgebra.Distributed.IterativeMethods;
	using MGroup.LinearAlgebra.Distributed.IterativeMethods.PCG;
	using MGroup.LinearAlgebra.Iterative.Termination.Iterations;

	public class FetiDPInterfaceProblemSolverFactoryPcg : IFetiDPInterfaceProblemSolverFactory
	{
		public int MaxIterations { get; set; } = 1000;

		public double ResidualTolerance { get; set; } = 1E-7;

		public bool ThrowExceptionIfNotConvergence { get; set; } = true;

		public bool UseObjectiveConvergenceCriterion { get; set; } = false;

		public IDistributedIterativeMethod BuildIterativeMethod(IPcgResidualConvergence convergenceCriterion)
		{
			var pcgBuilder = new PcgAlgorithm.Builder();
			pcgBuilder.ResidualTolerance = this.ResidualTolerance;
			pcgBuilder.MaxIterationsProvider = new FixedMaxIterationsProvider(this.MaxIterations);
			pcgBuilder.Convergence = convergenceCriterion;
			pcgBuilder.ThrowExceptionIfNotConvergence = this.ThrowExceptionIfNotConvergence;
			return pcgBuilder.Build();
		}
	}
}
