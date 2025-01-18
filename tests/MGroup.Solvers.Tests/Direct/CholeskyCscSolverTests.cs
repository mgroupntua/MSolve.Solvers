namespace MGroup.Solvers.Tests.Direct
{
	using System;
	using System.Collections.Generic;
	using System.Linq;
	using System.Text;
	using System.Threading.Tasks;

	using MGroup.Constitutive.Structural;
	using MGroup.LinearAlgebra.Implementations;
	using MGroup.LinearAlgebra.Implementations.NativeWin64;
	using MGroup.NumericalAnalyzers;
	using MGroup.Solvers.Direct;
	using MGroup.Solvers.DofOrdering;
	using MGroup.Solvers.DofOrdering.Reordering;
	using MGroup.Solvers.Tests.Benchmarks;
	using Xunit;

	public class CholeskyCscSolverTests
	{
		[Theory]
		[MemberData(nameof(TestSettings.ProvidersToTest), MemberType = typeof(TestSettings))]
		public static void TestWithoutReordering(IImplementationProvider provider)
		{
			CantileverBeam benchmark = new CantileverBeam.Builder().BuildWithQuad4Elements(200, 10);

			var solverFactory = new CholeskyCscSolver.Factory(provider);
			solverFactory.DofOrderer = new DofOrderer(new NodeMajorDofOrderingStrategy(), new NullReordering());
			var algebraicModel = solverFactory.BuildAlgebraicModel(benchmark.Model);
			using (CholeskyCscSolver solver = solverFactory.BuildSolver(algebraicModel))
			{
				benchmark.RunAnalysisAndCheck(algebraicModel, solver);
			}
		}

		[Theory]
		[MemberData(nameof(TestSettings.ProvidersToTest), MemberType = typeof(TestSettings))]
		public static void TestWithAmdReordering(IImplementationProvider provider)
		{
			CantileverBeam benchmark = new CantileverBeam.Builder().BuildWithQuad4Elements(200, 10);

			var solverFactory = new CholeskyCscSolver.Factory(provider);
			solverFactory.DofOrderer = new DofOrderer(new NodeMajorDofOrderingStrategy(), new AmdReordering(provider)); // default
			var algebraicModel = solverFactory.BuildAlgebraicModel(benchmark.Model);
			using (CholeskyCscSolver solver = solverFactory.BuildSolver(algebraicModel))
			{
				benchmark.RunAnalysisAndCheck(algebraicModel, solver);
			}
		}

		public static void SuiteSparseMemoryConsumptionDebugging()
		{
			for (int rep = 0; rep < 10; ++rep)
			{
				var benchmarkBuilder = new CantileverBeam.Builder();
				//benchmarkBuilder.Length = 5.0;
				CantileverBeam benchmark = benchmarkBuilder.BuildWithQuad4Elements(2000, 100);

				// Solver
				var provider = new NativeWin64ImplementationProvider();
				var solverFactory = new CholeskyCscSolver.Factory(provider);
				var algebraicModel = solverFactory.BuildAlgebraicModel(benchmark.Model);
				using (CholeskyCscSolver solver = solverFactory.BuildSolver(algebraicModel))
				{
					// Structural problem provider
					var problem = new ProblemStructural(benchmark.Model, algebraicModel);

					// Linear static analysis
					var childAnalyzer = new LinearAnalyzer(algebraicModel, solver, problem);
					var parentAnalyzer = new StaticAnalyzer(algebraicModel, problem, childAnalyzer);

					// Run the analysis
					parentAnalyzer.Initialize();
					parentAnalyzer.Solve();
				}
			}
		}
	}
}
