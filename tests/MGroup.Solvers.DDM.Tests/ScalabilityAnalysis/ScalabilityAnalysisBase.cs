namespace MGroup.Solvers.DDM.Tests.ScalabilityAnalysis
{
	using MGroup.Constitutive.Structural;
	using MGroup.Environments;
	using MGroup.LinearAlgebra.Implementations;
	using MGroup.MSolve.Discretization.Entities;
	using MGroup.MSolve.Solution;
	using MGroup.MSolve.Solution.AlgebraicModel;
	using MGroup.NumericalAnalyzers;

	public abstract class ScalabilityAnalysisBase
	{
		public IModelBuilder ModelBuilder { get; set; }

		public bool EnableNativeDlls { get; set; }

		public double IterativeResidualTolerance { get; set; } = 1E-6;

		public int NumSolverIterations { get; set; } = -1;

		//public int InterfaceProblemSize { get; set; } = -1;

		//public int CoarseProblemSize { get; set; } = -1;

		public void Clear()
		{
			NumSolverIterations = -1;
			//InterfaceProblemSize = -1;
			//CoarseProblemSize = -1;
		}

		public void RunParametricConstNumSubdomains(
			IComputeEnvironment environment, IImplementationProvider laProviderForSolver, string outputDirectory)
		{
			string path = outputDirectory + "results_const_subdomains.txt";
			environment.DoGlobalOperation(() => RenewFile(path));

			(List<int[]> numElements, int[] numSubdomains) = ModelBuilder.GetParametricConfigConstNumSubdomains();
			for (int i = 0; i < numElements.Count; i++)
			{
				Clear();
				ModelBuilder.NumElementsPerAxis = numElements[i];
				ModelBuilder.NumSubdomainsPerAxis = numSubdomains;
				RunSingleAnalysis(environment, laProviderForSolver);
				environment.DoGlobalOperation(() => PrintAnalysisData(path));
			}
		}

		public void RunParametricConstNumElements(
			IComputeEnvironment environment, IImplementationProvider laProviderForSolver, string outputDirectory)
		{
			string path = outputDirectory + "results_const_elements.txt";
			environment.DoGlobalOperation(() => RenewFile(path));

			(int[] numElements, List<int[]> numSubdomains) = ModelBuilder.GetParametricConfigConstNumElements();
			for (int i = 0; i < numSubdomains.Count; i++)
			{
				Clear();
				ModelBuilder.NumElementsPerAxis = numElements;
				ModelBuilder.NumSubdomainsPerAxis = numSubdomains[i];
				RunSingleAnalysis(environment, laProviderForSolver);
				environment.DoGlobalOperation(() => PrintAnalysisData(path));
			}
		}

		public void RunParametricConstSubdomainPerElementSize(
			IComputeEnvironment environment, IImplementationProvider laProviderForSolver, string outputDirectory)
		{

			string path = outputDirectory + "results_const_subdomain_per_element_size.txt";
			environment.DoGlobalOperation(() => RenewFile(path));

			(List<int[]> numElements, List<int[]> numSubdomains) = ModelBuilder.GetParametricConfigConstSubdomainPerElementSize();
			for (int i = 0; i < numSubdomains.Count; i++)
			{
				Clear();
				ModelBuilder.NumElementsPerAxis = numElements[i];
				ModelBuilder.NumSubdomainsPerAxis = numSubdomains[i];
				RunSingleAnalysis(environment, laProviderForSolver);
				environment.DoGlobalOperation(() => PrintAnalysisData(path));
			}
		}

		public void RunSingleAnalysis(IComputeEnvironment environment, IImplementationProvider laProviderForSolver)
		{
			(IModel model, ComputeNodeTopology nodeTopology) = ModelBuilder.CreateMultiSubdomainModel();
			model.ConnectDataStructures();
			(ISolver solver, IAlgebraicModel algebraicModel)
				= CreateSolver(environment, laProviderForSolver, model, nodeTopology);

			// Structural problem provider
			var provider = new ProblemStructural(model, algebraicModel);

			// Linear static analysis
			var childAnalyzer = new LinearAnalyzer(algebraicModel, solver, provider);
			var parentAnalyzer = new StaticAnalyzer(algebraicModel, provider, childAnalyzer);

			// Run the analysis
			parentAnalyzer.Initialize();
			parentAnalyzer.Solve();

			this.NumSolverIterations = solver.Logger.GetNumIterationsOfIterativeAlgorithm(0);
			//this.InterfaceProblemSize = solver.Logger.GetNumDofs(0, "Global boundary dofs");
			//this.CoarseProblemSize = solver.Logger.GetNumDofs(0, "Global corner dofs");
		}

		public abstract (ISolver solver, IAlgebraicModel algebraicModel) CreateSolver(IComputeEnvironment environment,
			IImplementationProvider laProviderForSolver, IModel model, ComputeNodeTopology nodeTopology);

		private void PrintAnalysisData(string path)
		{
			using (var writer = new StreamWriter(path, true))
			{
				writer.WriteLine("###############################################");
				writer.WriteLine(DateTime.Now);

				writer.Write("Number of elements per axis: ");
				for (int j = 0; j < ModelBuilder.NumElementsPerAxis.Length; j++)
				{
					writer.Write(ModelBuilder.NumElementsPerAxis[j] + " ");
				}
				writer.WriteLine();

				writer.Write("Number of subdomains per axis: ");
				for (int j = 0; j < ModelBuilder.NumSubdomainsPerAxis.Length; j++)
				{
					writer.Write(ModelBuilder.NumSubdomainsPerAxis[j] + " ");
				}
				writer.WriteLine();

				writer.Write("Subdomains size / element size: ");
				for (int j = 0; j < ModelBuilder.SubdomainSizePerElementSize.Length; j++)
				{
					writer.Write(ModelBuilder.SubdomainSizePerElementSize[j] + " ");
				}
				writer.WriteLine();

				//writer.WriteLine($"Interface problem size: {InterfaceProblemSize}");
				//writer.WriteLine($"Coarse problem size: {CoarseProblemSize}");
				writer.WriteLine($"Number of solver iterations: {NumSolverIterations}");
				writer.WriteLine();
			}
		}

		private void RenewFile(string path)
		{
			if (File.Exists(path))
			{
				File.Delete(path);
			}
		}
	}
}
