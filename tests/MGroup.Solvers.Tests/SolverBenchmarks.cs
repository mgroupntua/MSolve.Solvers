namespace MGroup.Solvers.Tests
{
	using MGroup.Constitutive.Structural;
	using MGroup.NumericalAnalyzers;
	using MGroup.Solvers.Direct;
	using MGroup.Solvers.Tests.Benchmarks;

	public class SolverBenchmarks
	{
		public static void SuiteSparseMemoryConsumptionDebugging()
		{
			for (int rep = 0; rep < 10; ++rep)
			{
				var benchmarkBuilder = new CantileverBeam.Builder();
				//benchmarkBuilder.Length = 5.0;
				CantileverBeam benchmark = benchmarkBuilder.BuildWithQuad4Elements(2000, 100);

				// Solver
				var solverFactory = new CholeskyCscSolver.Factory();
				var algebraicModel = solverFactory.BuildAlgebraicModel(benchmark.Model);
				using (CholeskyCscSolver solver = solverFactory.BuildSolver(algebraicModel))
				{
					// Structural problem provider
					var provider = new ProblemStructural(benchmark.Model, algebraicModel);

					// Linear static analysis
					var childAnalyzer = new LinearAnalyzer(algebraicModel, solver, provider);
					var parentAnalyzer = new StaticAnalyzer(algebraicModel, provider, childAnalyzer);

					// Run the analysis
					parentAnalyzer.Initialize();
					parentAnalyzer.Solve();
				}
			}
		}
	}
}
