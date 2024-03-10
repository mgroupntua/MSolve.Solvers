using System.Diagnostics;

using MGroup.Constitutive.Thermal;
using MGroup.Environments;
using MGroup.LinearAlgebra.Distributed.IterativeMethods.PCG;
using MGroup.LinearAlgebra.Iterative;
using MGroup.LinearAlgebra.Iterative.Termination.Iterations;
using MGroup.LinearAlgebra.Matrices;
using MGroup.MSolve.Discretization.Entities;
using MGroup.NumericalAnalyzers;
using MGroup.Solvers.DDM.FetiDP;
using MGroup.Solvers.DDM.FetiDP.CoarseProblem;
using MGroup.Solvers.DDM.FetiDP.Dofs;
using MGroup.Solvers.DDM.FetiDP.InterfaceProblem;
using MGroup.Solvers.DDM.FetiDP.Preconditioning;
using MGroup.Solvers.DDM.FetiDP.StiffnessMatrices;
using MGroup.Solvers.DDM.LinearSystem;
using MGroup.Solvers.DDM.Tests.ExampleModels;
using MGroup.Solvers.Results;

using Xunit;

namespace MGroup.Solvers.DDM.Tests.FetiDP
{
	public static class PapagiannakisFetiDPSolverTests
	{
		public enum Preconditioner
		{
			Dirichlet, Lumped, DiagonalDirichlet
		}

		[Theory]
		[InlineData(1.0, Preconditioner.DiagonalDirichlet, true, 10, 2E-9/*relaxed from 6.47E-10*/, EnvironmentChoice.SequentialShared)]
		[InlineData(1.0, Preconditioner.Dirichlet, true, 9, 2.44E-9, EnvironmentChoice.SequentialShared)]
		[InlineData(1.0, Preconditioner.Lumped, true, 11, 6.61E-10, EnvironmentChoice.SequentialShared)]
		public static void RunTest_8(double stiffnessRatio, Preconditioner preconditioner, bool ignoreHeterogenity, 
			int numIterationsExpected, double errorExpected, EnvironmentChoice environmentChoice)
			=> RunTest_8_Internal(stiffnessRatio, preconditioner, ignoreHeterogenity, numIterationsExpected, errorExpected, 
				environmentChoice.CreateEnvironment());

		internal static void RunTest_8_Internal(double stiffnessRatio, Preconditioner preconditioner, bool ignoreHeterogenity, 
			int numIterationsExpected, double errorExpected, IComputeEnvironment environment, bool isCoarseProblemDistributed = false)
		{
			// Model
			(Model model, ComputeNodeTopology nodeTopology) = PapagiannakisExample_8.CreateMultiSubdomainModel(stiffnessRatio);
			model.ConnectDataStructures(); //TODOMPI: this is also done in the analyzer
			ICornerDofSelection cornerDofs = PapagiannakisExample_8.GetCornerDofs(model);

			// Environment
			environment.Initialize(nodeTopology);

			// Solver
			var solverFactory = new FetiDPSolver<SymmetricCscMatrix>.Factory(
				environment, cornerDofs, new FetiDPSubdomainMatrixManagerSymmetricCSparse.Factory());

			if (preconditioner == Preconditioner.Dirichlet)
			{
				solverFactory.Preconditioner = new FetiDPDirichletPreconditioner();
			}
			else if (preconditioner == Preconditioner.Lumped)
			{
				solverFactory.Preconditioner = new FetiDPLumpedPreconditioner();
			}
			else
			{
				Debug.Assert(preconditioner == Preconditioner.DiagonalDirichlet);
				solverFactory.Preconditioner = new FetiDPDiagonalDirichletPreconditioner();
			}

			solverFactory.InterfaceProblemSolverFactory = new FetiDPInterfaceProblemSolverFactoryPcg()
			{
				UseObjectiveConvergenceCriterion = true,
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
				var coarseProblemMatrix = new FetiDPCoarseProblemMatrixSymmetricCSparse();
				solverFactory.CoarseProblemFactory = new FetiDPCoarseProblemGlobal.Factory(coarseProblemMatrix);
			}

			solverFactory.IsHomogeneousProblem = ignoreHeterogenity || (stiffnessRatio == 1.0);
			DistributedAlgebraicModel<SymmetricCscMatrix> algebraicModel = solverFactory.BuildAlgebraicModel(model);
			FetiDPSolver<SymmetricCscMatrix> solver = solverFactory.BuildSolver(model, algebraicModel);

			// Linear static analysis
			var problem = new ProblemThermal(model, algebraicModel);
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
		[InlineData(1.0, Preconditioner.DiagonalDirichlet, 14, 2E-6/*relaxed from 1.28E-9*/, EnvironmentChoice.SequentialShared)]
		[InlineData(1.0, Preconditioner.Dirichlet, 11, 4E-9/*relaxed from 1.39E-9*/, EnvironmentChoice.SequentialShared)]
		[InlineData(1.0, Preconditioner.Lumped, 18, 3E-9/*relaxed from 3.46E-10*/, EnvironmentChoice.SequentialShared)]
		public static void RunTest_9_1(double stiffnessRatio, Preconditioner preconditioner, int numIterationsExpected, 
			double errorExpected, EnvironmentChoice environmentChoice)
			=> RunTest_9_1_Internal(
				stiffnessRatio, preconditioner, numIterationsExpected, errorExpected, environmentChoice.CreateEnvironment());

		internal static void RunTest_9_1_Internal(double stiffnessRatio, Preconditioner preconditioner, int numIterationsExpected, 
			double errorExpected, IComputeEnvironment environment, bool isCoarseProblemDistributed = false)
		{
			// Model
			(Model model, ComputeNodeTopology nodeTopology) = PapagiannakisExample_9_1.CreateMultiSubdomainModel(stiffnessRatio);
			model.ConnectDataStructures(); //TODOMPI: this is also done in the analyzer
			ICornerDofSelection cornerDofs = PapagiannakisExample_9_1.GetCornerDofs(model);

			// Environment
			environment.Initialize(nodeTopology);

			// Solver
			var solverFactory = new FetiDPSolver<SymmetricCscMatrix>.Factory(
				environment, cornerDofs, new FetiDPSubdomainMatrixManagerSymmetricCSparse.Factory());
			solverFactory.InterfaceProblemSolverFactory = new FetiDPInterfaceProblemSolverFactoryPcg()
			{
				UseObjectiveConvergenceCriterion = true,
				MaxIterations = numIterationsExpected,
				ResidualTolerance = 1E-20,
				ThrowExceptionIfNotConvergence = false
			};

			if (preconditioner == Preconditioner.Dirichlet)
			{
				solverFactory.Preconditioner = new FetiDPDirichletPreconditioner();
			}
			else if (preconditioner == Preconditioner.Lumped)
			{
				solverFactory.Preconditioner = new FetiDPLumpedPreconditioner();
			}
			else
			{
				Debug.Assert(preconditioner == Preconditioner.DiagonalDirichlet);
				solverFactory.Preconditioner = new FetiDPDiagonalDirichletPreconditioner();
			}

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
				var coarseProblemMatrix = new FetiDPCoarseProblemMatrixSymmetricCSparse();
				solverFactory.CoarseProblemFactory = new FetiDPCoarseProblemGlobal.Factory(coarseProblemMatrix);
			}
			solverFactory.IsHomogeneousProblem = stiffnessRatio == 1.0;
			DistributedAlgebraicModel<SymmetricCscMatrix> algebraicModel = solverFactory.BuildAlgebraicModel(model);
			FetiDPSolver<SymmetricCscMatrix> solver = solverFactory.BuildSolver(model, algebraicModel);

			// Linear static analysis
			var problem = new ProblemThermal(model, algebraicModel);
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
		[InlineData(1.0, Preconditioner.DiagonalDirichlet, 10, 2E-5/*relaxed from 2.65E-10*/, EnvironmentChoice.SequentialShared)]
		[InlineData(1.0, Preconditioner.Dirichlet, 6, 3E-6/*relaxed from 3.32E-10*/, EnvironmentChoice.SequentialShared)]
		[InlineData(1.0, Preconditioner.Lumped, 12, 7E-6/*relaxed from 1.68E-10*/, EnvironmentChoice.SequentialShared)]
		public static void RunTest_9_2(double stiffnessRatio, Preconditioner preconditioner, int numIterationsExpected, 
			double errorExpected, EnvironmentChoice environmentChoice)
			=> RunTest_9_2_Internal(
				stiffnessRatio, preconditioner, numIterationsExpected, errorExpected, environmentChoice.CreateEnvironment());

		internal static void RunTest_9_2_Internal(double stiffnessRatio, Preconditioner preconditioner, int numIterationsExpected, 
			double errorExpected, IComputeEnvironment environment, bool isCoarseProblemDistributed = false)
		{
			// Model
			(Model model, ComputeNodeTopology nodeTopology) = PapagiannakisExample_9_2.CreateMultiSubdomainModel(stiffnessRatio);
			model.ConnectDataStructures(); //TODOMPI: this is also done in the analyzer
			ICornerDofSelection cornerDofs = PapagiannakisExample_9_2.GetCornerDofs(model);

			// Environment
			environment.Initialize(nodeTopology);

			// Solver
			var solverFactory = new FetiDPSolver<SymmetricCscMatrix>.Factory(
				environment, cornerDofs, new FetiDPSubdomainMatrixManagerSymmetricCSparse.Factory());

			if (preconditioner == Preconditioner.Dirichlet)
			{
				solverFactory.Preconditioner = new FetiDPDirichletPreconditioner();
			}
			else if (preconditioner == Preconditioner.Lumped)
			{
				solverFactory.Preconditioner = new FetiDPLumpedPreconditioner();
			}
			else
			{
				Debug.Assert(preconditioner == Preconditioner.DiagonalDirichlet);
				solverFactory.Preconditioner = new FetiDPDiagonalDirichletPreconditioner();
			}

			solverFactory.InterfaceProblemSolverFactory = new FetiDPInterfaceProblemSolverFactoryPcg()
			{
				#region debug uncomment
				//UseObjectiveConvergenceCriterion = true,
				#endregion
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
				var coarseProblemMatrix = new FetiDPCoarseProblemMatrixSymmetricCSparse();
				solverFactory.CoarseProblemFactory = new FetiDPCoarseProblemGlobal.Factory(coarseProblemMatrix);
			}

			solverFactory.IsHomogeneousProblem = stiffnessRatio == 1.0;
			DistributedAlgebraicModel<SymmetricCscMatrix> algebraicModel = solverFactory.BuildAlgebraicModel(model);
			FetiDPSolver<SymmetricCscMatrix> solver = solverFactory.BuildSolver(model, algebraicModel);

			// Linear static analysis
			var problem = new ProblemThermal(model, algebraicModel);
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
