//TODO: add performance logging for solvers and gather all these in the same method.
namespace MGroup.Solvers.Tests
{
	using MGroup.Constitutive.Structural;
	using MGroup.LinearAlgebra;
	using MGroup.LinearAlgebra.Matrices;
	using MGroup.MSolve.DataStructures;
	using MGroup.MSolve.Solution.AlgebraicModel;
	using MGroup.NumericalAnalyzers;
	using MGroup.Solvers.Direct;
	using MGroup.Solvers.DofOrdering;
	using MGroup.Solvers.DofOrdering.Reordering;
	using MGroup.Solvers.Iterative;
	using Xunit;

	public static class SingleSubdomainTests
	{
		[SkippableFact]
		internal static void TestDenseSolver()
		{
			Skip.IfNot(TestSettings.TestMkl, TestSettings.MessageWhenSkippingMKL);

			// Dense solver is too slow for a ~17.000 dof linear system, without MKL
			TestSettings.RunMultiproviderTest(LinearAlgebraProviderChoice.MKL, delegate ()
			{
				CantileverBeam benchmark = BuildCantileverBenchmark();

				var solverFactory = new DenseMatrixSolver.Factory();
				var algebraicModel = solverFactory.BuildAlgebraicModel(benchmark.Model);
				DenseMatrixSolver solver = solverFactory.BuildSolver(algebraicModel);
				//solverFactory.DofOrderer = new DofOrderer(new NodeMajorDofOrderingStrategy(), new NullReordering()); // default

				RunAnalysisAndCheck(benchmark, algebraicModel, solver);
			});
		}

		[Fact]
		internal static void TestPcgJacobiSolver()
		{
			CantileverBeam benchmark = BuildCantileverBenchmark();

			//LibrarySettings.LinearAlgebraProviders = LinearAlgebraProviderChoice.MKL;
			var solverFactory = new PcgSolver.Factory();
			//solverFactory.DofOrderer = new DofOrderer(new NodeMajorDofOrderingStrategy(), new NullReordering()); // default
			var algebraicModel = solverFactory.BuildAlgebraicModel(benchmark.Model);
			PcgSolver solver = solverFactory.BuildSolver(algebraicModel);

			RunAnalysisAndCheck(benchmark, algebraicModel, solver);
		}

		[Fact]
		internal static void TestPcgJacobiSolverWithAmdReordering()
		{
			CantileverBeam benchmark = BuildCantileverBenchmark();

			var solverFactory = new PcgSolver.Factory();
			solverFactory.DofOrderer = new DofOrderer(
				new NodeMajorDofOrderingStrategy(), AmdReordering.CreateWithCSparseAmd());
			var algebraicModel = solverFactory.BuildAlgebraicModel(benchmark.Model);
			PcgSolver solver = solverFactory.BuildSolver(algebraicModel);

			RunAnalysisAndCheck(benchmark, algebraicModel, solver);
		}

		[Fact]
		internal static void TestSkylineSolver()
		{
			CantileverBeam benchmark = BuildCantileverBenchmark();

			var solverFactory = new SkylineSolver.Factory();
			//solverFactory.DofOrderer = new DofOrderer(new NodeMajorDofOrderingStrategy(), new NullReordering()); // default
			var algebraicModel = solverFactory.BuildAlgebraicModel(benchmark.Model);
			SkylineSolver solver = solverFactory.BuildSolver(algebraicModel);
			solver.LinearSystem = algebraicModel.LinearSystem;

			RunAnalysisAndCheck(benchmark, algebraicModel, solver);
		}

		[Fact]
		internal static void TestSkylineSolverWithAmdReordering()
		{
			CantileverBeam benchmark = BuildCantileverBenchmark();

			var solverFactory = new SkylineSolver.Factory();
			solverFactory.DofOrderer = new DofOrderer(
				new NodeMajorDofOrderingStrategy(), AmdReordering.CreateWithCSparseAmd());
			var algebraicModel = solverFactory.BuildAlgebraicModel(benchmark.Model);
			SkylineSolver solver = solverFactory.BuildSolver(algebraicModel);
			RunAnalysisAndCheck(benchmark, algebraicModel, solver);

		}

		[SkippableFact]
		internal static void TestSuiteSparseSolver()
		{
			Skip.IfNot(TestSettings.TestSuiteSparse, TestSettings.MessageWhenSkippingSuiteSparse);

			CantileverBeam benchmark = BuildCantileverBenchmark();

			var solverFactory = new SuiteSparseSolver.Factory();
			solverFactory.DofOrderer = new DofOrderer(new NodeMajorDofOrderingStrategy(), new NullReordering());
			var algebraicModel = solverFactory.BuildAlgebraicModel(benchmark.Model);
			using (SuiteSparseSolver solver = solverFactory.BuildSolver(algebraicModel))
			{
				RunAnalysisAndCheck(benchmark, algebraicModel, solver);
			}
		}

		[SkippableFact]
		internal static void TestSuiteSparseSolverWithAmdReordering()
		{
			Skip.IfNot(TestSettings.TestSuiteSparse, TestSettings.MessageWhenSkippingSuiteSparse);

			CantileverBeam benchmark = BuildCantileverBenchmark();

			var solverFactory = new SuiteSparseSolver.Factory();
			//solverFactory.DofOrderer = new DofOrderer(
			//    new NodeMajorDofOrderingStrategy(), AmdReordering.CreateWithSuiteSparseAmd()); // default
			var algebraicModel = solverFactory.BuildAlgebraicModel(benchmark.Model);
			using (SuiteSparseSolver solver = solverFactory.BuildSolver(algebraicModel))
			{
				RunAnalysisAndCheck(benchmark, algebraicModel, solver);
			}
		}

		private static CantileverBeam BuildCantileverBenchmark()
		{
			var benchmarkBuilder = new CantileverBeam.Builder();
			//benchmarkBuilder.Length = 5.0;
			return benchmarkBuilder.BuildWithQuad4Elements(200, 10);
		}

		private static void RunAnalysisAndCheck<TMatrix>(CantileverBeam benchmark, IAlgebraicModel algebraicModel,
			SingleSubdomainSolverBase<TMatrix> solver) 
			where TMatrix : class, IMatrix
		{
			// Structural problem provider
			var provider = new ProblemStructural(benchmark.Model, algebraicModel);

			// Linear static analysis
			var childAnalyzer = new LinearAnalyzer(algebraicModel, solver, provider);
			var parentAnalyzer = new StaticAnalyzer(algebraicModel, provider, childAnalyzer);

			// Run the analysis
			parentAnalyzer.Initialize();
			parentAnalyzer.Solve();

			// Check output
			double endDeflectionExpected = benchmark.CalculateEndDeflectionWithEulerBeamTheory();
			double endDeflectionComputed =
				benchmark.CalculateAverageEndDeflectionFromSolution(solver.LinearSystem.Solution, algebraicModel);
			var comparer = new ValueComparer(1E-2);
			Assert.True(comparer.AreEqual(endDeflectionExpected, endDeflectionComputed));
		}
	}
}
