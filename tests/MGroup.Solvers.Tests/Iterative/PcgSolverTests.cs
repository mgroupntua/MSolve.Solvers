namespace MGroup.Solvers.Tests.Iterative
{
	using System;
	using System.Collections.Generic;
	using System.Linq;
	using System.Text;
	using System.Threading.Tasks;

	using MGroup.LinearAlgebra.Implementations.Managed;
	using MGroup.Solvers.DofOrdering;
	using MGroup.Solvers.DofOrdering.Reordering;
	using MGroup.Solvers.Iterative;
	using MGroup.Solvers.Tests.Benchmarks;
	using Xunit;

	public class PcgSolverTests
	{
		[Fact]
		public static void TestWithJacobiPreconditioner()
		{
			CantileverBeam benchmark = new CantileverBeam.Builder().BuildWithQuad4Elements(200, 10);

			var solverFactory = new PcgSolver.Factory();
			//solverFactory.DofOrderer = new DofOrderer(new NodeMajorDofOrderingStrategy(), new NullReordering()); // default
			var algebraicModel = solverFactory.BuildAlgebraicModel(benchmark.Model);
			PcgSolver solver = solverFactory.BuildSolver(algebraicModel);

			benchmark.RunAnalysisAndCheck(algebraicModel, solver);
		}

		[Fact]
		public static void TestWithJacobiPreconditionerAmdReordering()
		{
			CantileverBeam benchmark = new CantileverBeam.Builder().BuildWithQuad4Elements(200, 10);

			var solverFactory = new PcgSolver.Factory();
			var amd = new AmdReordering(new ManagedSequentialImplementationProvider());
			solverFactory.DofOrderer = new DofOrderer(new NodeMajorDofOrderingStrategy(), amd);
			var algebraicModel = solverFactory.BuildAlgebraicModel(benchmark.Model);
			PcgSolver solver = solverFactory.BuildSolver(algebraicModel);

			benchmark.RunAnalysisAndCheck(algebraicModel, solver);
		}
	}
}
