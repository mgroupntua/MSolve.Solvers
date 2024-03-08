using MGroup.LinearAlgebra.Distributed.IterativeMethods;
using MGroup.LinearAlgebra.Distributed.IterativeMethods.PCG;

namespace MGroup.Solvers.DDM.FetiDP.InterfaceProblem
{
	public interface IFetiDPInterfaceProblemSolverFactory
	{
		int MaxIterations { get; set; }

		double ResidualTolerance { get; set; }

		bool ThrowExceptionIfNotConvergence { get; set; }

		bool UseObjectiveConvergenceCriterion { get; set; }

		IDistributedIterativeMethod BuildIterativeMethod(IPcgResidualConvergence convergenceCriterion);
	}
}
