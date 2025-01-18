namespace MGroup.Solvers.Tests.Direct
{
	using System;
	using System.Collections.Generic;
	using System.Linq;
	using System.Text;

	using MGroup.LinearAlgebra.Implementations.Managed;
	using MGroup.Solvers.Direct;
	using MGroup.Solvers.DofOrdering;
	using MGroup.Solvers.DofOrdering.Reordering;
	using MGroup.Solvers.Tests.Benchmarks;
	using Xunit;

	public class LUCscSolverTests
	{
		[Fact]
		public static void TestWithoutReordering()
		{
			CantileverBeam benchmark = new CantileverBeam.Builder().BuildWithQuad4Elements(200, 10);

			var solverFactory = new LUCscSolver.Factory();
			solverFactory.DofOrderer = new DofOrderer(new NodeMajorDofOrderingStrategy(), new NullReordering()); // default
			var algebraicModel = solverFactory.BuildAlgebraicModel(benchmark.Model);
			LUCscSolver solver = solverFactory.BuildSolver(algebraicModel);
			solver.LinearSystem = algebraicModel.LinearSystem;

			benchmark.RunAnalysisAndCheck(algebraicModel, solver);
		}

		[Fact]
		internal static void TestWithAmdReordering()
		{
			CantileverBeam benchmark = new CantileverBeam.Builder().BuildWithQuad4Elements(200, 10);

			var solverFactory = new LUCscSolver.Factory();
			var amd = new AmdReordering(new ManagedSequentialImplementationProvider());
			solverFactory.DofOrderer = new DofOrderer(new NodeMajorDofOrderingStrategy(), amd);
			var algebraicModel = solverFactory.BuildAlgebraicModel(benchmark.Model);
			LUCscSolver solver = solverFactory.BuildSolver(algebraicModel);

			benchmark.RunAnalysisAndCheck(algebraicModel, solver);
		}
	}
}
