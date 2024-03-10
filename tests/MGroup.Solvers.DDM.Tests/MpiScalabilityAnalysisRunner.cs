using MGroup.Environments.Mpi;
using MGroup.Solvers.DDM.Tests.ScalabilityAnalysis;

namespace MGroup.Solvers.DDM.Tests
{
	public class MpiScalabilityAnalysisRunner
	{
		public static void RunScalabilityAnalysesWith4Processes()
		{
			using (var mpiEnvironment = new MpiEnvironment(new MasterSlavesGlobalOperationStrategy()))
			{
				MpiDebugUtilities.AssistDebuggerAttachment();

				MpiDebugUtilities.DoSerially(MPI.Communicator.world,
					() => Console.WriteLine(
						$"Process {MPI.Communicator.world.Rank}: Now running ScalabilityAnalysisPFetiDP.RunFullScalabilityAnalysisCantilever2DInternal"));
				ScalabilityAnalysisPFetiDP.RunFullScalabilityAnalysisCantilever2DInternal(mpiEnvironment);

				MpiDebugUtilities.DoSerially(MPI.Communicator.world,
					() => Console.WriteLine(
						$"Process {MPI.Communicator.world.Rank}: Now running ScalabilityAnalysisPFetiDP.RunFullScalabilityAnalysisRve2DInternal"));
				ScalabilityAnalysisPFetiDP.RunFullScalabilityAnalysisRve2DInternal(mpiEnvironment);

				MpiDebugUtilities.DoSerially(MPI.Communicator.world,
					() => Console.WriteLine(
						$"Process {MPI.Communicator.world.Rank}: Now running ScalabilityAnalysisPFetiDP.RunFullScalabilityAnalysisCantilever3DInternal"));
				ScalabilityAnalysisPFetiDP.RunFullScalabilityAnalysisCantilever3DInternal(mpiEnvironment);

				MpiDebugUtilities.DoSerially(MPI.Communicator.world,
					() => Console.WriteLine(
						$"Process {MPI.Communicator.world.Rank}: Now running ScalabilityAnalysisPFetiDP.RunFullScalabilityAnalysisRve3DInternal"));
				ScalabilityAnalysisPFetiDP.RunFullScalabilityAnalysisRve3DInternal(mpiEnvironment);



				MpiDebugUtilities.DoSerially(MPI.Communicator.world,
					() => Console.WriteLine($"Process {MPI.Communicator.world.Rank}: All tests passed"));
			}
		}
	}
}
