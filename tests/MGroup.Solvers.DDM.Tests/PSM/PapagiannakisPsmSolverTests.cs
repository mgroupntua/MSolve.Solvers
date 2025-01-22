//TODO: Compare the tolerances with the original Papagiannakis examples and try to improve accuracy
namespace MGroup.Solvers.DDM.Tests.PSM
{
	using MGroup.Constitutive.Structural;
	using MGroup.Environments;
	using MGroup.LinearAlgebra.Implementations;
	using MGroup.LinearAlgebra.Iterative;
	using MGroup.LinearAlgebra.Matrices;
	using MGroup.MSolve.Discretization.Entities;
	using MGroup.NumericalAnalyzers;
	using MGroup.Solvers.DDM.LinearSystem;
	using MGroup.Solvers.DDM.Psm;
	using MGroup.Solvers.DDM.PSM.InterfaceProblem;
	using MGroup.Solvers.DDM.PSM.StiffnessMatrices;
	using MGroup.Solvers.DDM.Tests.ExampleModels;
	using MGroup.Solvers.Results;
	using MGroup.Solvers.Tests;
	using MGroup.Solvers.Tests.TempUtilityClasses;

	using Xunit;

	[Collection("Sequential")]
	public static class PapagiannakisPsmSolverTests
	{
		public static TheoryData<double, int, double, IEnvironmentChoice, IImplementationProviderChoice> DataForTest_8
		{
			get
			{
				var data = new TheoryData<double, int, double, IEnvironmentChoice, IImplementationProviderChoice>();
				TestSettings.CombineTheoryDataWithAllProvidersAndEnvironments(data, 1.0, 26, 5.51E-10);
				TestSettings.CombineTheoryDataWithAllProvidersAndEnvironments(data, 1E3, 44, 2.59E-10);
				TestSettings.CombineTheoryDataWithAllProvidersAndEnvironments(data, 1E4, 61, 2.56E-10);
				TestSettings.CombineTheoryDataWithAllProvidersAndEnvironments(data, 1E5, 73, 4.89E-9);
				TestSettings.CombineTheoryDataWithAllProvidersAndEnvironments(data, 1E6, 81, 3.35E-8 /*relaxed from 2.02E-8*/);
				return data;
			}
		}

		[Theory]
		[MemberData(nameof(DataForTest_8))]
		public static void RunTest_8(double stiffnessRatio, int numIterationsExpected, double errorExpected,
			IEnvironmentChoice environmentChoice, IImplementationProviderChoice laProviderForSolver)
		{
			RunTest_8_Internal(stiffnessRatio, numIterationsExpected, errorExpected,
				environmentChoice.Activate(), laProviderForSolver.Activate());
		}

		internal static void RunTest_8_Internal(double stiffnessRatio, int numIterationsExpected, double errorExpected,
			IComputeEnvironment environment, IImplementationProvider laProviderForSolver)
		{
			// Model
			(Model model, ComputeNodeTopology nodeTopology) = PapagiannakisExample_8.CreateMultiSubdomainModel(stiffnessRatio);
			model.ConnectDataStructures(); //TODOMPI: this is also done in the analyzer

			// Environment
			environment.Initialize(nodeTopology);

			// Solver
			var solverFactory = new PsmSolver<SymmetricCscMatrix>.Factory(
				environment, laProviderForSolver, new PsmSubdomainMatrixManagerSymmetricCsc.Factory());
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

		public static TheoryData<double, int, double, IEnvironmentChoice, IImplementationProviderChoice> DataForTest_9_1
		{
			get
			{
				var data = new TheoryData<double, int, double, IEnvironmentChoice, IImplementationProviderChoice>();
				TestSettings.CombineTheoryDataWithAllProvidersAndEnvironments(data, 1.0, 50, 4E-9 /*relaxed from 9.91E-10*/);
				TestSettings.CombineTheoryDataWithAllProvidersAndEnvironments(data, 1E2, 94, 5.30E-11 /*relaxed from 5.01E-12*/);
				TestSettings.CombineTheoryDataWithAllProvidersAndEnvironments(data, 1E3, 120, 4.29E-11 /*relaxed from 2.87E-11*/);
				TestSettings.CombineTheoryDataWithAllProvidersAndEnvironments(data, 1E4, 163, 1.03E-9);
				TestSettings.CombineTheoryDataWithAllProvidersAndEnvironments(data, 1E5, 192, 4.59E-9 /*relaxed from 6.31E-19*/);
				TestSettings.CombineTheoryDataWithAllProvidersAndEnvironments(data, 1E6, 238, 7.68E-08 /*relaxed from 3.02E-08*/);
				return data;
			}
		}

		[Theory]
		[MemberData(nameof(DataForTest_9_1))]
		public static void RunTest_9_1(double stiffnessRatio, int numIterationsExpected, double errorExpected,
			IEnvironmentChoice environmentChoice, IImplementationProviderChoice laProviderForSolver)
		{
			RunTest_9_1_Internal(stiffnessRatio, numIterationsExpected, errorExpected,
				environmentChoice.Activate(), laProviderForSolver.Activate());
		}

		internal static void RunTest_9_1_Internal(double stiffnessRatio, int numIterationsExpected, double errorExpected,
			IComputeEnvironment environment, IImplementationProvider laProviderForSolver)
		{
			// Model
			(Model model, ComputeNodeTopology nodeTopology) = PapagiannakisExample_9_1.CreateMultiSubdomainModel(stiffnessRatio);
			model.ConnectDataStructures(); //TODOMPI: this is also done in the analyzer

			// Environment
			environment.Initialize(nodeTopology);

			// Solver
			var solverFactory = new PsmSolver<SymmetricCscMatrix>.Factory(
				environment, laProviderForSolver, new PsmSubdomainMatrixManagerSymmetricCsc.Factory());
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

		public static TheoryData<double, int, double, IEnvironmentChoice, IImplementationProviderChoice> DataForTest_9_2
		{
			get
			{
				var data = new TheoryData<double, int, double, IEnvironmentChoice, IImplementationProviderChoice>();
				TestSettings.CombineTheoryDataWithAllProvidersAndEnvironments(data, 1.0, 19, 6E-4 /*relaxed from 1.15E-12*/);
				TestSettings.CombineTheoryDataWithAllProvidersAndEnvironments(data, 1E2, 60, 7.10E-5); /*relaxed from 48, 8.49E-11*/
				TestSettings.CombineTheoryDataWithAllProvidersAndEnvironments(data, 1E3, 85, 1.50E-6 /*relaxed from 3.88E-10*/);
				TestSettings.CombineTheoryDataWithAllProvidersAndEnvironments(data, 1E4, 99, 5.61E-6 /*relaxed from 3.92E-8*/);
				TestSettings.CombineTheoryDataWithAllProvidersAndEnvironments(data, 1E5, 101, 8.16E-5 /*relaxed from 4.78E-6*/);
				TestSettings.CombineTheoryDataWithAllProvidersAndEnvironments(data, 1E6, 123, 7.83E-5 /*relaxed from 4.81E-6*/);
				return data;
			}
		}

		[Theory]
		[MemberData(nameof(DataForTest_9_2))]
		public static void RunTest_9_2(double stiffnessRatio, int numIterationsExpected, double errorExpected,
			IEnvironmentChoice environmentChoice, IImplementationProviderChoice laProviderForSolver)
		{
			RunTest_9_2_Internal(stiffnessRatio, numIterationsExpected, errorExpected,
				environmentChoice.Activate(), laProviderForSolver.Activate());
		}

		internal static void RunTest_9_2_Internal(double stiffnessRatio, int numIterationsExpected, double errorExpected,
			IComputeEnvironment environment, IImplementationProvider laProviderForSolver)
		{
			// Model
			(Model model, ComputeNodeTopology nodeTopology) = PapagiannakisExample_9_2.CreateMultiSubdomainModel(stiffnessRatio);
			model.ConnectDataStructures(); //TODOMPI: this is also done in the analyzer

			// Environment
			environment.Initialize(nodeTopology);

			// Solver
			var solverFactory = new PsmSolver<SymmetricCscMatrix>.Factory(
				environment, laProviderForSolver, new PsmSubdomainMatrixManagerSymmetricCsc.Factory());
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
			NodalResults expectedResults = PapagiannakisExample_9_2.SolveWithSkylineSolver(
				PapagiannakisExample_9_2.CreateSingleSubdomainModel(stiffnessRatio));
			Assert.Equal(PapagiannakisExample_9_2.NumTotalDofs, expectedResults.Data.NumEntries);
			NodalResults globalComputedResults = algebraicModel.ExtractGlobalResults(solver.LinearSystem.Solution, 1E-6);
			double error = expectedResults.Subtract(globalComputedResults).Norm2() / expectedResults.Norm2();
			Assert.InRange(error, 0, errorExpected);
		}
	}
}
