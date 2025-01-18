namespace MGroup.Solvers.Tests.Direct
{
	using System;
	using System.Collections.Generic;
	using System.Linq;
	using System.Text;
	using System.Threading.Tasks;

	using MGroup.LinearAlgebra.Implementations.NativeWin64;
	using MGroup.Solvers.Direct;
	using MGroup.Solvers.DofOrdering;
	using MGroup.Solvers.DofOrdering.Reordering;
	using MGroup.Solvers.Tests.Benchmarks;

	using Xunit;

	public class DenseSolverTests
	{
		[SkippableFact]
		internal static void TestDenseSolver()
		{
			// Dense solver is too slow for a ~17.000 dof linear system, without MKL
			Skip.IfNot(TestSettings.LibsToTest.Win64IntelMkl, TestSettings.SkipMessage);
			LinearAlgebra.LibrarySettings.GlobalProvider = new NativeWin64ImplementationProvider();

			CantileverBeam benchmark = new CantileverBeam.Builder().BuildWithQuad4Elements(200, 10);

			var solverFactory = new DenseMatrixSolver.Factory();
			var algebraicModel = solverFactory.BuildAlgebraicModel(benchmark.Model);
			DenseMatrixSolver solver = solverFactory.BuildSolver(algebraicModel);
			solverFactory.DofOrderer = new DofOrderer(new NodeMajorDofOrderingStrategy(), new NullReordering()); // default

			benchmark.RunAnalysisAndCheck(algebraicModel, solver);
		}
	}
}
