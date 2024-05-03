namespace MGroup.Solvers.DDM.PSM.InterfaceProblem
{
	using MGroup.LinearAlgebra.Distributed.IterativeMethods;
	using MGroup.LinearAlgebra.Distributed.IterativeMethods.PCG;

	public interface IPsmInterfaceProblemSolverFactory
	{
		int MaxIterations { get; set; }

		double ResidualTolerance { get; set; }

		bool ThrowExceptionIfNotConvergence { get; set; }

		bool UseObjectiveConvergenceCriterion { get; set; }

		IDistributedIterativeMethod BuildIterativeMethod(IPcgResidualConvergence convergenceCriterion);
	}
}
