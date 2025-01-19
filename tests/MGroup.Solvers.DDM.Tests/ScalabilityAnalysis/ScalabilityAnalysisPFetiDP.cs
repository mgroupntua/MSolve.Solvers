namespace MGroup.Solvers.DDM.Tests.ScalabilityAnalysis
{
	using MGroup.Environments;
	using MGroup.Environments.Mpi;
	using MGroup.LinearAlgebra.Distributed.IterativeMethods.PCG;
	using MGroup.LinearAlgebra.Implementations;
	using MGroup.LinearAlgebra.Implementations.NativeWin64;
	using MGroup.LinearAlgebra.Iterative.Termination.Iterations;
	using MGroup.LinearAlgebra.Matrices;
	using MGroup.MSolve.Discretization.Entities;
	using MGroup.MSolve.Solution;
	using MGroup.MSolve.Solution.AlgebraicModel;
	using MGroup.Solvers.DDM.FetiDP.CoarseProblem;
	using MGroup.Solvers.DDM.FetiDP.Dofs;
	using MGroup.Solvers.DDM.FetiDP.StiffnessMatrices;
	using MGroup.Solvers.DDM.LinearSystem;
	using MGroup.Solvers.DDM.PFetiDP;
	using MGroup.Solvers.DDM.Psm;
	using MGroup.Solvers.DDM.PSM.InterfaceProblem;
	using MGroup.Solvers.DDM.PSM.StiffnessMatrices;

	public class ScalabilityAnalysisPFetiDP : ScalabilityAnalysisBase
	{
		private const string workingDirectory = @"C:\Users\Serafeim\Desktop\PFETIDP\scalability\";

		//[Theory]
		//[InlineData(EnvironmentChoice.SequentialSharedEnvironment)]
		//[InlineData(EnvironmentChoice.TplSharedEnvironment)]
		public static void RunFullScalabilityAnalysisCantilever2D(EnvironmentChoice environmentChoice)
			=> RunFullScalabilityAnalysisCantilever2DInternal(environmentChoice.CreateEnvironment());

		internal static void RunFullScalabilityAnalysisCantilever2DInternal(IComputeEnvironment environment)
		{
			var laProviderForSolver = new NativeWin64ImplementationProvider();

			string outputDirectory = workingDirectory + @"\cantilever2D\";
			var scalabilityAnalysis = new ScalabilityAnalysisPFetiDP();
			scalabilityAnalysis.ModelBuilder = new CantilevelBeam2D();
			scalabilityAnalysis.EnableNativeDlls = true;
			scalabilityAnalysis.IterativeResidualTolerance = 1E-6;

			scalabilityAnalysis.RunParametricConstNumSubdomains(environment, laProviderForSolver, outputDirectory);
			//scalabilityAnalysis.RunParametricConstNumElements(environment, laProviderForSolver, outputDirectory);
			//scalabilityAnalysis.RunParametricConstSubdomainPerElementSize(environment, laProviderForSolver, outputDirectory);
		}

		//[Theory]
		//[InlineData(EnvironmentChoice.SequentialSharedEnvironment)]
		//[InlineData(EnvironmentChoice.TplSharedEnvironment)]
		public static void RunFullScalabilityAnalysisCantilever3D(EnvironmentChoice environmentChoice)
			=> RunFullScalabilityAnalysisCantilever3DInternal(environmentChoice.CreateEnvironment());

		internal static void RunFullScalabilityAnalysisCantilever3DInternal(IComputeEnvironment environment)
		{
			var laProviderForSolver = new NativeWin64ImplementationProvider();

			string outputDirectory = workingDirectory + @"\cantilever3D\";
			var scalabilityAnalysis = new ScalabilityAnalysisPFetiDP();
			scalabilityAnalysis.ModelBuilder = new CantilevelBeam3D();
			scalabilityAnalysis.EnableNativeDlls = true;
			scalabilityAnalysis.IterativeResidualTolerance = 1E-6;

			scalabilityAnalysis.RunParametricConstNumSubdomains(environment, laProviderForSolver, outputDirectory);
			//scalabilityAnalysis.RunParametricConstNumElements(environment, laProviderForSolver, outputDirectory);
			//scalabilityAnalysis.RunParametricConstSubdomainPerElementSize(environment, laProviderForSolver, outputDirectory);
		}

		//[Theory]
		//[InlineData(EnvironmentChoice.SequentialSharedEnvironment)]
		//[InlineData(EnvironmentChoice.TplSharedEnvironment)]
		public static void RunFullScalabilityAnalysisRve2D(EnvironmentChoice environmentChoice)
			=> RunFullScalabilityAnalysisRve2DInternal(environmentChoice.CreateEnvironment());

		internal static void RunFullScalabilityAnalysisRve2DInternal(IComputeEnvironment environment)
		{
			var laProviderForSolver = new NativeWin64ImplementationProvider();

			string outputDirectory = workingDirectory + @"\rve2D\";
			var scalabilityAnalysis = new ScalabilityAnalysisPFetiDP();
			scalabilityAnalysis.ModelBuilder = new Rve2D();
			scalabilityAnalysis.EnableNativeDlls = true;
			scalabilityAnalysis.IterativeResidualTolerance = 1E-6;

			scalabilityAnalysis.RunParametricConstNumSubdomains(environment, laProviderForSolver, outputDirectory);
			//scalabilityAnalysis.RunParametricConstNumElements(environment, laProviderForSolver, outputDirectory);
			//scalabilityAnalysis.RunParametricConstSubdomainPerElementSize(environment, laProviderForSolver, outputDirectory);
		}

		//[Theory]
		//[InlineData(EnvironmentChoice.SequentialSharedEnvironment)]
		//[InlineData(EnvironmentChoice.TplSharedEnvironment)]
		public static void RunFullScalabilityAnalysisRve3D(EnvironmentChoice environmentChoice)
			=> RunFullScalabilityAnalysisRve3DInternal(environmentChoice.CreateEnvironment());

		internal static void RunFullScalabilityAnalysisRve3DInternal(IComputeEnvironment environment)
		{
			var laProviderForSolver = new NativeWin64ImplementationProvider();

			string outputDirectory = workingDirectory + @"\rve3D\";
			var scalabilityAnalysis = new ScalabilityAnalysisPFetiDP();
			scalabilityAnalysis.ModelBuilder = new Rve3D();
			scalabilityAnalysis.EnableNativeDlls = true;
			scalabilityAnalysis.IterativeResidualTolerance = 1E-6;

			scalabilityAnalysis.RunParametricConstNumSubdomains(environment, laProviderForSolver, outputDirectory);
			//scalabilityAnalysis.RunParametricConstNumElements(environment, laProviderForSolver, outputDirectory);
			//scalabilityAnalysis.RunParametricConstSubdomainPerElementSize(environment, laProviderForSolver, outputDirectory);
		}

		public override (ISolver solver, IAlgebraicModel algebraicModel) CreateSolver(IComputeEnvironment environment,
			IImplementationProvider laProviderForSolver, IModel model, ComputeNodeTopology nodeTopology)
		{
			ICornerDofSelection cornerDofs = ModelBuilder.GetCornerDofs(model);
			environment.Initialize(nodeTopology);

			// Specify the format of matrices
			var psmMatricesFactory = new PsmSubdomainMatrixManagerSymmetricCsc.Factory();
			var fetiDPMatricesFactory = new FetiDPSubdomainMatrixManagerSymmetricCsc.Factory(true);
			var coarseProblemMatrix = new FetiDPCoarseProblemMatrixSymmetricCsc(laProviderForSolver);

			var solverFactory = new PFetiDPSolver<SymmetricCscMatrix>.Factory(
				environment, laProviderForSolver, psmMatricesFactory, cornerDofs, fetiDPMatricesFactory);

			solverFactory.InterfaceProblemSolverFactory = new PsmInterfaceProblemSolverFactoryPcg()
			{
				MaxIterations = 200,
				ResidualTolerance = 1E-7
			};

			//TODOMPI: This was written before global coarse problem was implemented for MPI. 
			//		I should probably use that instead of distributed coarse problem.
			if (environment is MpiEnvironment) 
			{
				var coarseProblemPcgBuilder = new PcgAlgorithm.Builder();
				coarseProblemPcgBuilder.MaxIterationsProvider = new FixedMaxIterationsProvider(200);
				coarseProblemPcgBuilder.ResidualTolerance = 1E-12;
				var coarseProblemFactory = new FetiDPCoarseProblemDistributed.Factory();
				coarseProblemFactory.CoarseProblemSolver = coarseProblemPcgBuilder.Build();
				solverFactory.CoarseProblemFactory = coarseProblemFactory;
			}
			else
			{
				solverFactory.CoarseProblemFactory = new FetiDPCoarseProblemGlobal.Factory(coarseProblemMatrix);
			}

			solverFactory.IsHomogeneousProblem = true;
			DistributedAlgebraicModel<SymmetricCscMatrix> algebraicModel = solverFactory.BuildAlgebraicModel(model);
			PsmSolver<SymmetricCscMatrix> solver = solverFactory.BuildSolver(model, algebraicModel);

			return (solver, algebraicModel);
		}
	}
}
