namespace MGroup.Solvers.DDM.Tests
{
	using MGroup.Environments;
	using MGroup.Environments.Mpi;
	using MGroup.LinearAlgebra.Implementations;
	using MGroup.Solvers.DDM.Tests.ScalabilityAnalysis;
	using MGroup.Solvers.Tests.TempUtilityClasses;
	using MGroup.Solvers.Tests;
	using MGroup.Solvers.DDM.Tests.PSM;

	public class MpiScalabilityAnalysisRunner
	{
		public static void RunScalabilityAnalysesWith4Processes()
		{
			using (var mpiEnvironment = new MpiEnvironment(new MasterSlavesGlobalOperationStrategy()))
			{
				MpiDebugUtilities.AssistDebuggerAttachment();

				RunTestInMpiEnvironmentWithAllProviders(mpiEnvironment,
					"ScalabilityAnalysisPFetiDP.RunFullScalabilityAnalysisCantilever2DInternal",
					ScalabilityAnalysisPFetiDP.RunFullScalabilityAnalysisCantilever2DInternal);

				RunTestInMpiEnvironmentWithAllProviders(mpiEnvironment,
					"ScalabilityAnalysisPFetiDP.RunFullScalabilityAnalysisRve2DInternal",
					ScalabilityAnalysisPFetiDP.RunFullScalabilityAnalysisRve2DInternal);

				RunTestInMpiEnvironmentWithAllProviders(mpiEnvironment,
					"ScalabilityAnalysisPFetiDP.RunFullScalabilityAnalysisCantilever3DInternal",
					ScalabilityAnalysisPFetiDP.RunFullScalabilityAnalysisCantilever3DInternal);

				RunTestInMpiEnvironmentWithAllProviders(mpiEnvironment,
					"ScalabilityAnalysisPFetiDP.RunFullScalabilityAnalysisRve3DInternal",
					ScalabilityAnalysisPFetiDP.RunFullScalabilityAnalysisRve3DInternal);

				MpiDebugUtilities.DoSerially(MPI.Communicator.world,
					() => Console.WriteLine($"Process {MPI.Communicator.world.Rank}: All tests passed"));
			}
		}

		private static void RunTestInMpiEnvironmentWithAllProviders(MpiEnvironment mpiEnvironment, string testName,
			Action<IComputeEnvironment, IImplementationProvider> test)
		{
			//TODO: get the name of the test by the delegate (probably using Reflection) instead of having it as an extra param
			foreach (IImplementationProviderChoice provider in TestSettings.ProviderChoicesToTest)
			{
				MpiDebugUtilities.DoSerially(MPI.Communicator.world,
					() => Console.WriteLine(
						$"Process {MPI.Communicator.world.Rank} on processor {MPI.Environment.ProcessorName}: " +
						$"Now running example {testName}, with provider {provider.ToString()}."));
				test(mpiEnvironment, provider.Activate());
			}
		}
	}
}
