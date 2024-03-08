using System;
using System.Collections.Generic;
using System.Text;
using MGroup.LinearAlgebra.Distributed.IterativeMethods;
using MGroup.LinearAlgebra.Distributed.IterativeMethods.PCG;

namespace MGroup.Solvers.DDM.PSM.InterfaceProblem
{
	public interface IPsmInterfaceProblemSolverFactory
	{
		int MaxIterations { get; set; }

		double ResidualTolerance { get; set; }

		bool ThrowExceptionIfNotConvergence { get; set; }

		bool UseObjectiveConvergenceCriterion { get; set; }

		IDistributedIterativeMethod BuildIterativeMethod(IPcgResidualConvergence convergenceCriterion);
	}
}
