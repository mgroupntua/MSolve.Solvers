//TODO: Allow client to register tests to be run, schedule them, run them with appropriate before and after messages, 
//      handle exceptions by adequately notifying user, but without crushing and omitting the next tests. 
//TODO: Hande the case of 1 machine running multiple tests that need different number of processes (make sure the max needed 
//      have been launched, deactivate the extras per test so that they do not mess with collectives). 
//      Handle the case of multiple machines running tests. Facilite user with command line parameters to MPI exec
namespace MGroup.Solvers.DDM.Tests
{
	using MGroup.Environments;
	using MGroup.Environments.Mpi;
	using MGroup.LinearAlgebra.Implementations;
	using MGroup.LinearAlgebra.Implementations.Managed;
	using MGroup.Solvers.DDM.Tests.PFetiDP;
	using MGroup.Solvers.DDM.Tests.PSM;
	using MGroup.Solvers.Tests;
	using MGroup.Solvers.Tests.TempUtilityClasses;

	using Microsoft.VisualStudio.TestPlatform.PlatformAbstractions.Interfaces;

	public static class MpiTestSuite
	{
		public static void RunTestsWith5Processes()
		{
			IMpiGlobalOperationStrategy globalOperationStrategy = new MasterSlavesGlobalOperationStrategy(0);

			// Process 4 will be used only for global operations
			//IMpiGlobalOperationStrategy globalOperationStrategy = new MasterSlavesGlobalOperationStrategy(4);

			//IMpiGlobalOperationStrategy globalOperationStrategy = new DemocraticGlobalOperationStrategy();
			using (var mpiEnvironment = new MpiEnvironment(globalOperationStrategy))
			{
				MpiDebugUtilities.AssistDebuggerAttachment();

				RunTestInMpiEnvironment(mpiEnvironment, "PsmInterfaceProblemDofsTests.TestForLine1DInternal",
					PsmInterfaceProblemDofsTests.TestForLine1DInternal);

				RunTestInMpiEnvironment(mpiEnvironment, "PsmInterfaceProblemDofsTests.TestForPlane2DInternal",
					PsmInterfaceProblemDofsTests.TestForPlane2DInternal);

				RunTestInMpiEnvironmentWithAllProviders(mpiEnvironment, "SimplePsmSolverTests.TestForLine1DInternal",
					SimplePsmSolverTests.TestForLine1DInternal);

				RunTestInMpiEnvironmentWithAllProviders(mpiEnvironment, "SimplePsmSolverTests.TestForPlane2DInternal",
					SimplePsmSolverTests.TestForPlane2DInternal);

				RunTestInMpiEnvironmentWithAllProviders(mpiEnvironment, "SimplePsmSolverTests.TestForBrick3DInternal",
					SimplePsmSolverTests.TestForBrick3DInternal);

				RunTestInMpiEnvironmentWithAllProviders(mpiEnvironment, 
					"SimplePFetiDPSolverTests.TestForPlane2DInternal with distributed coarse problem",
					SimplePFetiDPSolverTests.TestForPlane2DInternal, true, false, false);

				RunTestInMpiEnvironmentWithAllProviders(mpiEnvironment,
					"SimplePFetiDPSolverTests.TestForBrick3DInternal with distributed coarse problem",
					SimplePFetiDPSolverTests.TestForBrick3DInternal, true, false, false);

				RunTestInMpiEnvironmentWithAllProviders(mpiEnvironment,
					"SimplePFetiDPSolverTests.TestForPlane2DInternal with global coarse problem",
					SimplePFetiDPSolverTests.TestForPlane2DInternal, false, false, false);

				RunTestInMpiEnvironmentWithAllProviders(mpiEnvironment,
					"SimplePFetiDPSolverTests.TestForBrick3DInternal with global coarse problem",
					SimplePFetiDPSolverTests.TestForBrick3DInternal, false, false, false);

				MpiDebugUtilities.DoSerially(MPI.Communicator.world,
					() => Console.WriteLine(
						$"Process {MPI.Communicator.world.Rank} on processor {MPI.Environment.ProcessorName}: " +
						$"All tests passed"));
			}
		}

		private static void RunTestInMpiEnvironment(MpiEnvironment mpiEnvironment, string testName,
			Action<IComputeEnvironment> test)
		{
			MpiDebugUtilities.DoSerially(MPI.Communicator.world,
					() => Console.WriteLine(
						$"Process {MPI.Communicator.world.Rank} on processor {MPI.Environment.ProcessorName}: " +
						$"Now running test {testName}."));
			test(mpiEnvironment);
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
						$"Now running test {testName}, with provider {provider.ToString()}."));
				test(mpiEnvironment, provider.Activate());
			}
		}

		private static void RunTestInMpiEnvironmentWithAllProviders<T1>(MpiEnvironment mpiEnvironment, string testName,
			Action<T1, IComputeEnvironment, IImplementationProvider> test, T1 param1)
		{
			//TODO: get the name of the test by the delegate (probably using Reflection) instead of having it as an extra param
			foreach (IImplementationProviderChoice provider in TestSettings.ProviderChoicesToTest)
			{
				MpiDebugUtilities.DoSerially(MPI.Communicator.world,
					() => Console.WriteLine(
						$"Process {MPI.Communicator.world.Rank} on processor {MPI.Environment.ProcessorName}: " +
						$"Now running test {testName}, with provider {provider.ToString()}."));
				test(param1, mpiEnvironment, provider.Activate());
			}
		}

		private static void RunTestInMpiEnvironmentWithAllProviders<T1, T2>(MpiEnvironment mpiEnvironment, string testName,
			Action<T1, T2, IComputeEnvironment, IImplementationProvider> test, T1 param1, T2 param2)
		{
			//TODO: get the name of the test by the delegate (probably using Reflection) instead of having it as an extra param
			foreach (IImplementationProviderChoice provider in TestSettings.ProviderChoicesToTest)
			{
				MpiDebugUtilities.DoSerially(MPI.Communicator.world,
					() => Console.WriteLine(
						$"Process {MPI.Communicator.world.Rank} on processor {MPI.Environment.ProcessorName}: " +
						$"Now running test {testName}, with provider {provider.ToString()}."));
				test(param1, param2, mpiEnvironment, provider.Activate());
			}
		}

		private static void RunTestInMpiEnvironmentWithAllProviders<T1, T2, T3>(MpiEnvironment mpiEnvironment, string testName,
			Action<T1, T2, T3, IComputeEnvironment, IImplementationProvider> test, T1 param1, T2 param2, T3 param3)
		{
			//TODO: get the name of the test by the delegate (probably using Reflection) instead of having it as an extra param
			foreach (IImplementationProviderChoice provider in TestSettings.ProviderChoicesToTest)
			{
				MpiDebugUtilities.DoSerially(MPI.Communicator.world,
					() => Console.WriteLine(
						$"Process {MPI.Communicator.world.Rank} on processor {MPI.Environment.ProcessorName}: " +
						$"Now running test {testName}, with provider {provider.ToString()}."));
				test(param1, param2, param3, mpiEnvironment, provider.Activate());
			}
		}
	}
}
