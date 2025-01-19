namespace MGroup.Solvers.DDM.Tests.PFetiDP
{
	using MGroup.Constitutive.Structural;
	using MGroup.Environments;
	using MGroup.LinearAlgebra.Distributed.IterativeMethods.PCG;
	using MGroup.LinearAlgebra.Implementations.Managed;
	using MGroup.LinearAlgebra.Iterative;
	using MGroup.LinearAlgebra.Iterative.Termination.Iterations;
	using MGroup.LinearAlgebra.Matrices;
	using MGroup.MSolve.Discretization.Entities;
	using MGroup.NumericalAnalyzers;
	using MGroup.Solvers.DDM.FetiDP.CoarseProblem;
	using MGroup.Solvers.DDM.FetiDP.Dofs;
	using MGroup.Solvers.DDM.FetiDP.StiffnessMatrices;
	using MGroup.Solvers.DDM.LinearSystem;
	using MGroup.Solvers.DDM.PFetiDP;
	using MGroup.Solvers.DDM.Psm;
	using MGroup.Solvers.DDM.PSM.InterfaceProblem;
	using MGroup.Solvers.DDM.PSM.StiffnessMatrices;
	using MGroup.Solvers.DDM.Tests.ExampleModels;
	using MGroup.Solvers.Results;

	using Xunit;

	[Collection("Sequential")]
	public static class PapagiannakisPFetiDPSolverTests
	{
		[Theory]
		[InlineData(1.0, true, 10, 1.53E-9, EnvironmentChoice.SequentialShared)]
		[InlineData(1E3, false, 11, 2.32E-10, EnvironmentChoice.SequentialShared)]
		[InlineData(1E3, true, 25, 2.86E-10, EnvironmentChoice.SequentialShared)]
		[InlineData(1E4, false, 11, 3.47E-10 /*relaxed from 1.73E-10*/, EnvironmentChoice.SequentialShared)]
		[InlineData(1E4, true, 33, 1.46E-9, EnvironmentChoice.SequentialShared)]
		[InlineData(1E5, false, 11, 4E-9 /*relaxed from  1.05E-9*/, EnvironmentChoice.SequentialShared)]
		[InlineData(1E5, true, 38, 3E-9 /*relaxed from 5.9E-10*/, EnvironmentChoice.SequentialShared)]
		[InlineData(1E6, false, 11, 2.00E-7, EnvironmentChoice.SequentialShared)]
		[InlineData(1E6, true, 53, 2.24E-7, EnvironmentChoice.SequentialShared)]
		public static void RunTest_8(double stiffnessRatio, bool ignoreHeterogenity, int numIterationsExpected, 
			double errorExpected, EnvironmentChoice environmentChoice)
			=> RunTest_8_Internal(stiffnessRatio, ignoreHeterogenity, numIterationsExpected, errorExpected, 
				environmentChoice.CreateEnvironment());

		internal static void RunTest_8_Internal(double stiffnessRatio, bool ignoreHeterogenity, int numIterationsExpected,
			double errorExpected, IComputeEnvironment environment, bool isCoarseProblemDistributed = false)
		{
			var laProviderForSolver = new ManagedSequentialImplementationProvider();

			// Model
			(Model model, ComputeNodeTopology nodeTopology) = PapagiannakisExample_8.CreateMultiSubdomainModel(stiffnessRatio);
			model.ConnectDataStructures(); //TODOMPI: this is also done in the analyzer
			ICornerDofSelection cornerDofs = PapagiannakisExample_8.GetCornerDofs(model);

			// Environment
			environment.Initialize(nodeTopology);

			// Solver
			var solverFactory = new PFetiDPSolver<SymmetricCscMatrix>.Factory(
				environment, laProviderForSolver, new PsmSubdomainMatrixManagerSymmetricCsc.Factory(),
				cornerDofs, new FetiDPSubdomainMatrixManagerSymmetricCsc.Factory(true));

			solverFactory.InterfaceProblemSolverFactory = new PsmInterfaceProblemSolverFactoryPcg()
			{
				// Papagiannakis specified these and reported the number of iterations and the error from direct solver.
				//MaxIterations = 1000,
				//ResidualTolerance = 1E-7,
				UseObjectiveConvergenceCriterion = true,

				// Instead I will set a tolerance that is impossible to reach, let PCG run for the same number of iterations
				// as Papagiannakis and compare the error from direct solver.
				MaxIterations = numIterationsExpected,
				ResidualTolerance = 1E-20,
				ThrowExceptionIfNotConvergence = false
			};

			if (isCoarseProblemDistributed)
			{
				var pcgBuilder = new PcgAlgorithm.Builder();
				pcgBuilder.MaxIterationsProvider = new FixedMaxIterationsProvider(200);
				pcgBuilder.ResidualTolerance = 1E-7;
				var coarseProblemFactory = new FetiDPCoarseProblemDistributed.Factory();
				coarseProblemFactory.CoarseProblemSolver = pcgBuilder.Build();
				solverFactory.CoarseProblemFactory = coarseProblemFactory;
			}
			else
			{
				var coarseProblemMatrix = new FetiDPCoarseProblemMatrixSymmetricCsc(laProviderForSolver);
				solverFactory.CoarseProblemFactory = new FetiDPCoarseProblemGlobal.Factory(coarseProblemMatrix);
			}

			solverFactory.IsHomogeneousProblem = ignoreHeterogenity || (stiffnessRatio == 1.0);
			DistributedAlgebraicModel<SymmetricCscMatrix> algebraicModel = solverFactory.BuildAlgebraicModel(model);
			PsmSolver<SymmetricCscMatrix> solver = solverFactory.BuildSolver(model, algebraicModel);

			// Linear static analysis
			var problem = new ProblemStructural(model, algebraicModel);
			var childAnalyzer = new LinearAnalyzer(algebraicModel, solver, problem);
			var parentAnalyzer = new StaticAnalyzer(algebraicModel, problem, childAnalyzer);

			// Run the analysis
			parentAnalyzer.Initialize();
			parentAnalyzer.Solve();

			// Check convergence
			IterativeStatistics stats = solver.InterfaceProblemSolutionStats;
			//Assert.InRange(stats.NumIterationsRequired, 1, numIterationsExpected); // Do not check this. It is guaranteed.

			// Check results
			NodalResults expectedResults = PapagiannakisExample_8.SolveWithSkylineSolver(
				PapagiannakisExample_8.CreateSingleSubdomainModel(stiffnessRatio));
			Assert.Equal(PapagiannakisExample_8.NumTotalDofs, expectedResults.Data.NumEntries);
			NodalResults globalComputedResults = algebraicModel.ExtractGlobalResults(solver.LinearSystem.Solution, 1E-6);
			double error = expectedResults.Subtract(globalComputedResults).Norm2() / expectedResults.Norm2();
			Assert.InRange(error, 0, errorExpected);
		}

		[Theory]
		[InlineData(1.0, 11, 2E-8 /*relaxed from 4.94E-9*/, EnvironmentChoice.SequentialShared)]
		[InlineData(1E2, 11, 7E-10 /*relaxed from 3.06E-10*/, EnvironmentChoice.SequentialShared)]
		[InlineData(1E3, 12, 5.22E-11 /*relaxed from 1.14E-11*/, EnvironmentChoice.SequentialShared)]
		[InlineData(1E4, 12, 9.92E-10, EnvironmentChoice.SequentialShared)]
		[InlineData(1E5, 12, 7.76E-9, EnvironmentChoice.SequentialShared)]
		[InlineData(1E6, 13, 1E-7 /*relaxed from 2.97E-8*/, EnvironmentChoice.SequentialShared)]
		public static void RunTest_9_1(
			double stiffnessRatio, int numIterationsExpected, double errorExpected, EnvironmentChoice environmentChoice)
			=> RunTest_9_1_Internal(
				stiffnessRatio, numIterationsExpected, errorExpected, environmentChoice.CreateEnvironment());

		internal static void RunTest_9_1_Internal(double stiffnessRatio, int numIterationsExpected, double errorExpected, 
			IComputeEnvironment environment, bool isCoarseProblemDistributed = false)
		{
			var laProviderForSolver = new ManagedSequentialImplementationProvider();

			// Model
			(Model model, ComputeNodeTopology nodeTopology) = PapagiannakisExample_9_1.CreateMultiSubdomainModel(stiffnessRatio);
			model.ConnectDataStructures(); //TODOMPI: this is also done in the analyzer
			ICornerDofSelection cornerDofs = PapagiannakisExample_9_1.GetCornerDofs(model);

			// Environment
			environment.Initialize(nodeTopology);

			// Solver
			var solverFactory = new PFetiDPSolver<SymmetricCscMatrix>.Factory(
				environment, laProviderForSolver, new PsmSubdomainMatrixManagerSymmetricCsc.Factory(),
				cornerDofs, new FetiDPSubdomainMatrixManagerSymmetricCsc.Factory(true));
			solverFactory.InterfaceProblemSolverFactory = new PsmInterfaceProblemSolverFactoryPcg()
			{
				// Papagiannakis specified these and reported the number of iterations and the error from direct solver.
				//MaxIterations = 1000,
				//ResidualTolerance = 1E-5,
				UseObjectiveConvergenceCriterion = true,

				// Instead I will set a tolerance that is impossible to reach, let PCG run for the same number of iterations
				// as Papagiannakis and compare the error from direct solver.
				MaxIterations = numIterationsExpected,
				ResidualTolerance = 1E-20,
				ThrowExceptionIfNotConvergence = false
			};

			if (isCoarseProblemDistributed)
			{
				var pcgBuilder = new PcgAlgorithm.Builder();
				pcgBuilder.MaxIterationsProvider = new FixedMaxIterationsProvider(200);
				pcgBuilder.ResidualTolerance = 1E-5;
				var coarseProblemFactory = new FetiDPCoarseProblemDistributed.Factory();
				coarseProblemFactory.CoarseProblemSolver = pcgBuilder.Build();
				solverFactory.CoarseProblemFactory = coarseProblemFactory;
			}
			else
			{
				var coarseProblemMatrix = new FetiDPCoarseProblemMatrixSymmetricCsc(laProviderForSolver);
				solverFactory.CoarseProblemFactory = new FetiDPCoarseProblemGlobal.Factory(coarseProblemMatrix);
			}
			solverFactory.IsHomogeneousProblem = stiffnessRatio == 1.0;
			DistributedAlgebraicModel<SymmetricCscMatrix> algebraicModel = solverFactory.BuildAlgebraicModel(model);
			PsmSolver<SymmetricCscMatrix> solver = solverFactory.BuildSolver(model, algebraicModel);

			// Linear static analysis
			var problem = new ProblemStructural(model, algebraicModel);
			var childAnalyzer = new LinearAnalyzer(algebraicModel, solver, problem);
			var parentAnalyzer = new StaticAnalyzer(algebraicModel, problem, childAnalyzer);

			// Run the analysis
			parentAnalyzer.Initialize();
			parentAnalyzer.Solve();

			// Check convergence
			IterativeStatistics stats = solver.InterfaceProblemSolutionStats;
			//Assert.InRange(stats.NumIterationsRequired, 1, numIterationsExpected); // Do not check this. It is guaranteed.

			// Check results
			NodalResults expectedResults = PapagiannakisExample_9_1.SolveWithSkylineSolver(
				PapagiannakisExample_9_1.CreateSingleSubdomainModel(stiffnessRatio));
			Assert.Equal(PapagiannakisExample_9_1.NumTotalDofs, expectedResults.Data.NumEntries);
			NodalResults globalComputedResults = algebraicModel.ExtractGlobalResults(solver.LinearSystem.Solution, 1E-6);
			double error = expectedResults.Subtract(globalComputedResults).Norm2() / expectedResults.Norm2();
			Assert.InRange(error, 0, errorExpected);
		}

		[Theory]
		[InlineData(1.0, 7, 1E-7 /*relaxed from 1.09E-10*/, EnvironmentChoice.SequentialShared)]
		[InlineData(1E2, 8, 2E-7 /*relaxed from 8.43E-10*/, EnvironmentChoice.SequentialShared)]
		[InlineData(1E3, 7, 6E-7 /*relaxed from 4.74E-8*/, EnvironmentChoice.SequentialShared)]
		[InlineData(1E4, 7, 5E-7 /*relaxed from 4.96E-8*/, EnvironmentChoice.SequentialShared)]
		[InlineData(1E5, 7, 5E-7 /*relaxed from 5.14E-8*/, EnvironmentChoice.SequentialShared)]
		[InlineData(1E6, 7, 5E-7 /*relaxed from 5.20E-8*/, EnvironmentChoice.SequentialShared)]
		public static void RunTest_9_2(
			double stiffnessRatio, int numIterationsExpected, double errorExpected, EnvironmentChoice environmentChoice)
			=> RunTest_9_2_Internal(
				stiffnessRatio, numIterationsExpected, errorExpected, environmentChoice.CreateEnvironment());

		internal static void RunTest_9_2_Internal(double stiffnessRatio, int numIterationsExpected, double errorExpected,
			IComputeEnvironment environment, bool isCoarseProblemDistributed = false)
		{
			var laProviderForSolver = new ManagedSequentialImplementationProvider();

			// Model
			(Model model, ComputeNodeTopology nodeTopology) = PapagiannakisExample_9_2.CreateMultiSubdomainModel(stiffnessRatio);
			model.ConnectDataStructures(); //TODOMPI: this is also done in the analyzer
			ICornerDofSelection cornerDofs = PapagiannakisExample_9_2.GetCornerDofs(model);

			// Environment
			environment.Initialize(nodeTopology);

			// Solver
			var solverFactory = new PFetiDPSolver<SymmetricCscMatrix>.Factory(
				environment, laProviderForSolver, new PsmSubdomainMatrixManagerSymmetricCsc.Factory(),
				cornerDofs, new FetiDPSubdomainMatrixManagerSymmetricCsc.Factory(true));

			solverFactory.InterfaceProblemSolverFactory = new PsmInterfaceProblemSolverFactoryPcg()
			{
				// Papagiannakis specified these and reported the number of iterations and the error from direct solver.
				//MaxIterations = 1000,
				//ResidualTolerance = 1E-7,
				UseObjectiveConvergenceCriterion = true,

				// Instead I will set a tolerance that is impossible to reach, let PCG run for the same number of iterations
				// as Papagiannakis and compare the error from direct solver.
				MaxIterations = numIterationsExpected,
				ResidualTolerance = 1E-20,
				ThrowExceptionIfNotConvergence = false
			};

			if (isCoarseProblemDistributed)
			{
				var pcgBuilder = new PcgAlgorithm.Builder();
				pcgBuilder.MaxIterationsProvider = new FixedMaxIterationsProvider(200);
				pcgBuilder.ResidualTolerance = 1E-7;
				var coarseProblemFactory = new FetiDPCoarseProblemDistributed.Factory();
				coarseProblemFactory.CoarseProblemSolver = pcgBuilder.Build();
				solverFactory.CoarseProblemFactory = coarseProblemFactory;
			}
			else
			{
				var coarseProblemMatrix = new FetiDPCoarseProblemMatrixSymmetricCsc(laProviderForSolver);
				solverFactory.CoarseProblemFactory = new FetiDPCoarseProblemGlobal.Factory(coarseProblemMatrix);
			}

			solverFactory.IsHomogeneousProblem = stiffnessRatio == 1.0;
			DistributedAlgebraicModel<SymmetricCscMatrix> algebraicModel = solverFactory.BuildAlgebraicModel(model);
			PsmSolver<SymmetricCscMatrix> solver = solverFactory.BuildSolver(model, algebraicModel);

			// Linear static analysis
			var problem = new ProblemStructural(model, algebraicModel);
			var childAnalyzer = new LinearAnalyzer(algebraicModel, solver, problem);
			var parentAnalyzer = new StaticAnalyzer(algebraicModel, problem, childAnalyzer);

			// Run the analysis
			parentAnalyzer.Initialize();
			parentAnalyzer.Solve();

			// Check convergence
			IterativeStatistics stats = solver.InterfaceProblemSolutionStats;
			int relaxedIterationsExpected = 2 + numIterationsExpected;
			//Assert.InRange(stats.NumIterationsRequired, 1, numIterationsExpected); // Do not check this. It is guaranteed.

			// Check results
			NodalResults expectedResults = PapagiannakisExample_9_2.SolveWithSkylineSolver(
				PapagiannakisExample_9_2.CreateSingleSubdomainModel(stiffnessRatio));
			Assert.Equal(PapagiannakisExample_9_2.NumTotalDofs, expectedResults.Data.NumEntries);
			NodalResults globalComputedResults = algebraicModel.ExtractGlobalResults(solver.LinearSystem.Solution, 1E-6);
			double error = expectedResults.Subtract(globalComputedResults).Norm2() / expectedResults.Norm2();
			Assert.InRange(error, 0, errorExpected);
		}
	}
}
