using MGroup.Constitutive.Thermal;
using MGroup.Environments;
using MGroup.LinearAlgebra.Distributed.IterativeMethods.PCG;
using MGroup.LinearAlgebra.Distributed.IterativeMethods.PCG.Reorthogonalization;
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
	public static class SimpleFetiDPSolverTests
	{
		[Theory]
		[InlineData(EnvironmentChoice.SequentialShared, false, false, false)]
		[InlineData(EnvironmentChoice.TplShared, false, false, false)]
		public static void TestForBrick3D(EnvironmentChoice env, bool coarseDistributed, bool coarseJacobi, bool coarseReortho)
			=> TestForBrick3DInternal(env.CreateEnvironment(), coarseDistributed, coarseJacobi, coarseReortho);

		internal static void TestForBrick3DInternal(IComputeEnvironment environment, bool isCoarseProblemDistributed,
			bool useCoarseJacobiPreconditioner, bool useReorthogonalizedPcg)
		{
			// Environment
			ComputeNodeTopology nodeTopology = Brick3DExample.CreateNodeTopology();
			environment.Initialize(nodeTopology);

			// Model
			IModel model = Brick3DExample.CreateMultiSubdomainModel();
			model.ConnectDataStructures(); //TODOMPI: this is also done in the analyzer
			ICornerDofSelection cornerDofs = Brick3DExample.GetCornerDofs(model);

			// Solver
			var solverFactory = new FetiDPSolver<SymmetricCscMatrix>.Factory(
				environment, cornerDofs, new FetiDPSubdomainMatrixManagerSymmetricCSparse.Factory());

			solverFactory.Preconditioner = new FetiDPDirichletPreconditioner();
			solverFactory.InterfaceProblemSolverFactory = new FetiDPInterfaceProblemSolverFactoryPcg()
			{
				MaxIterations = 200,
				ResidualTolerance = 1E-10
			};

			if (isCoarseProblemDistributed)
			{
				var coarseProblemFactory = new FetiDPCoarseProblemDistributed.Factory();
				if (useReorthogonalizedPcg)
				{
					var coarseProblemPcgBuilder = new ReorthogonalizedPcg.Builder();
					coarseProblemPcgBuilder.MaxIterationsProvider = new FixedMaxIterationsProvider(200);
					coarseProblemPcgBuilder.ResidualTolerance = 2E-12;
					coarseProblemPcgBuilder.DirectionVectorsRetention = new FixedDirectionVectorsRetention(40, true);
					coarseProblemFactory.CoarseProblemSolver = coarseProblemPcgBuilder.Build();
				}
				else
				{
					var coarseProblemPcgBuilder = new PcgAlgorithm.Builder();
					coarseProblemPcgBuilder.MaxIterationsProvider = new FixedMaxIterationsProvider(200);
					coarseProblemPcgBuilder.ResidualTolerance = 2E-12;
					coarseProblemFactory.CoarseProblemSolver = coarseProblemPcgBuilder.Build();
				}

				coarseProblemFactory.UseJacobiPreconditioner = useCoarseJacobiPreconditioner;
				solverFactory.CoarseProblemFactory = coarseProblemFactory;
			}
			else 
			{
				var coarseProblemMatrix = new FetiDPCoarseProblemMatrixSymmetricCSparse();
				solverFactory.CoarseProblemFactory = new FetiDPCoarseProblemGlobal.Factory(coarseProblemMatrix);
			}

			DistributedAlgebraicModel<SymmetricCscMatrix> algebraicModel = solverFactory.BuildAlgebraicModel(model);
			FetiDPSolver<SymmetricCscMatrix> solver = solverFactory.BuildSolver(model, algebraicModel);

			// Linear static analysis
			var problem = new ProblemThermal(model, algebraicModel);
			var childAnalyzer = new LinearAnalyzer(algebraicModel, solver, problem);
			var parentAnalyzer = new StaticAnalyzer(algebraicModel, problem, childAnalyzer);

			// Run the analysis
			parentAnalyzer.Initialize();
			parentAnalyzer.Solve();

			// Check results
			NodalResults expectedResults = Brick3DExample.GetExpectedNodalValues(problem.ActiveDofs);
			double tolerance = 1E-7;
			environment.DoPerNode(subdomainID =>
			{
				NodalResults computedResults = algebraicModel.ExtractAllResults(subdomainID, solver.LinearSystem.Solution);
				Assert.True(expectedResults.IsSuperSetOf(computedResults, tolerance, out string msg), msg);
			});

			//Debug.WriteLine($"Num PCG iterations = {solver.PcgStats.NumIterationsRequired}," +
			//    $" final residual norm ratio = {solver.PcgStats.ResidualNormRatioEstimation}");

			// Check convergence
			int pcgIterationsExpected = 30;
			double pcgResidualNormRatioExpected = 6E-11;
			IterativeStatistics stats = solver.InterfaceProblemSolutionStats;
			Assert.Equal(pcgIterationsExpected, stats.NumIterationsRequired);
			Assert.InRange(stats.ResidualNormRatioEstimation, 0, pcgResidualNormRatioExpected);
		}

		[Theory]
		[InlineData(EnvironmentChoice.SequentialShared, false, false, false)]
		//[InlineData(EnvironmentChoice.SequentialShared, true, false, false)]
		//[InlineData(EnvironmentChoice.SequentialShared, true, true, false)]
		//[InlineData(EnvironmentChoice.SequentialShared, true, false, true)]
		//[InlineData(EnvironmentChoice.SequentialShared, true, true, true)]
		[InlineData(EnvironmentChoice.TplShared, false, false, false)]
		//[InlineData(EnvironmentChoice.TplShared, true, false, false)]
		//[InlineData(EnvironmentChoice.TplShared, true, true, false)]
		//[InlineData(EnvironmentChoice.TplShared, true, false, true)]
		//[InlineData(EnvironmentChoice.TplShared, true, true, true)]
		public static void TestForPlane2D(EnvironmentChoice env, bool coarseDistributed, bool coarseJacobi, bool coarseReortho)
			=> TestForPlane2DInternal(env.CreateEnvironment(), coarseDistributed, coarseJacobi, coarseReortho);

		internal static void TestForPlane2DInternal(IComputeEnvironment environment, bool isCoarseProblemDistributed,
			bool useCoarseJacobiPreconditioner, bool useReorthogonalizedPcg)
		{
			// Environment
			ComputeNodeTopology nodeTopology = Plane2DExample.CreateNodeTopology();
			environment.Initialize(nodeTopology);

			// Model
			IModel model = Plane2DExample.CreateMultiSubdomainModel();
			model.ConnectDataStructures(); //TODOMPI: this is also done in the analyzer
			ICornerDofSelection cornerDofs = Plane2DExample.GetCornerDofs(model);

			// Solver
			var solverFactory = new FetiDPSolver<SymmetricCscMatrix>.Factory(
				environment, cornerDofs, new FetiDPSubdomainMatrixManagerSymmetricCSparse.Factory());

			solverFactory.Preconditioner = new FetiDPDirichletPreconditioner();
			solverFactory.InterfaceProblemSolverFactory = new FetiDPInterfaceProblemSolverFactoryPcg()
			{
				MaxIterations = 200,
				ResidualTolerance = 1E-10
			};

			if (isCoarseProblemDistributed)
			{
				var coarseProblemFactory = new FetiDPCoarseProblemDistributed.Factory();
				if (useReorthogonalizedPcg)
				{
					var coarseProblemPcgBuilder = new ReorthogonalizedPcg.Builder();
					coarseProblemPcgBuilder.DirectionVectorsRetention = new FixedDirectionVectorsRetention(30, true);
					coarseProblemPcgBuilder.MaxIterationsProvider = new FixedMaxIterationsProvider(200);
					coarseProblemPcgBuilder.ResidualTolerance = 2E-12;
					coarseProblemFactory.CoarseProblemSolver = coarseProblemPcgBuilder.Build();
				}
				else
				{
					var coarseProblemPcgBuilder = new PcgAlgorithm.Builder();
					coarseProblemPcgBuilder.MaxIterationsProvider = new FixedMaxIterationsProvider(200);
					coarseProblemPcgBuilder.ResidualTolerance = 2E-12;
					coarseProblemFactory.CoarseProblemSolver = coarseProblemPcgBuilder.Build();
				}
				
				coarseProblemFactory.UseJacobiPreconditioner = useCoarseJacobiPreconditioner;
				solverFactory.CoarseProblemFactory = coarseProblemFactory;
			}
			else 
			{
				var coarseProblemMatrix = new FetiDPCoarseProblemMatrixSymmetricCSparse();
				solverFactory.CoarseProblemFactory = new FetiDPCoarseProblemGlobal.Factory(coarseProblemMatrix);
			}

			DistributedAlgebraicModel<SymmetricCscMatrix> algebraicModel = solverFactory.BuildAlgebraicModel(model);
			FetiDPSolver<SymmetricCscMatrix> solver = solverFactory.BuildSolver(model, algebraicModel);

			// Linear static analysis
			var problem = new ProblemThermal(model, algebraicModel);
			var childAnalyzer = new LinearAnalyzer(algebraicModel, solver, problem);
			var parentAnalyzer = new StaticAnalyzer(algebraicModel, problem, childAnalyzer);

			// Run the analysis
			parentAnalyzer.Initialize();
			parentAnalyzer.Solve();

			// Check results
			NodalResults expectedResults = Plane2DExample.GetExpectedNodalValues(problem.ActiveDofs);
			double tolerance = 1E-7;
			environment.DoPerNode(subdomainID =>
			{
				NodalResults computedResults = algebraicModel.ExtractAllResults(subdomainID, solver.LinearSystem.Solution);
				Assert.True(expectedResults.IsSuperSetOf(computedResults, tolerance, out string msg), msg);
			});

			//Debug.WriteLine($"Num PCG iterations = {solver.PcgStats.NumIterationsRequired}," +
			//    $" final residual norm ratio = {solver.PcgStats.ResidualNormRatioEstimation}");

			// Check convergence
			int pcgIterationsExpected = 14;
			double pcgResidualNormRatioExpected = 4E-11;
			IterativeStatistics stats = solver.InterfaceProblemSolutionStats;
			Assert.Equal(pcgIterationsExpected, stats.NumIterationsRequired);
			Assert.InRange(stats.ResidualNormRatioEstimation, 0, pcgResidualNormRatioExpected);
		}
	}
}
