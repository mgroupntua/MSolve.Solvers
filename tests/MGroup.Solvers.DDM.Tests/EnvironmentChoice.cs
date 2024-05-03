namespace MGroup.Solvers.DDM.Tests
{
	using System;

	using MGroup.Environments;
	using MGroup.Environments.Mpi;

	public enum EnvironmentChoice
	{
		SequentialShared, TplShared, Mpi
	}

	public static class EnvironmentChoiceExtensions
	{
		public static IComputeEnvironment CreateEnvironment(this EnvironmentChoice environmentChoice)
		{
			if (environmentChoice == EnvironmentChoice.SequentialShared)
			{
				return new SequentialSharedEnvironment();
			}
			else if (environmentChoice == EnvironmentChoice.TplShared)
			{
				return new TplSharedEnvironment();
			}
			else if (environmentChoice == EnvironmentChoice.Mpi)
			{
				return new MpiEnvironment(new MasterSlavesGlobalOperationStrategy());
			}
			else
			{
				throw new NotImplementedException();
			}
		}
	}
}
