namespace MGroup.Solvers.Tests.TempUtilityClasses
{
	using System;
	using System.Collections.Generic;
	using System.Linq;
	using System.Text;
	using System.Threading.Tasks;

	using MGroup.Environments;
	using MGroup.Environments.Mpi;

	using Xunit.Abstractions;

	public class MpiEnvironmentChoice : IEnvironmentChoice
	{
		private IMpiGlobalOperationStrategy globalOperationStrategy; //TODO: This needs its own serialization

		public MpiEnvironmentChoice()
		{
			this.globalOperationStrategy = new MasterSlavesGlobalOperationStrategy();
		}

		//public MpiEnvironmentChoice(IMpiGlobalOperationStrategy globalOperationStrategy)
		//{
		//	this.globalOperationStrategy = globalOperationStrategy;
		//}

		public IComputeEnvironment Activate() => new MpiEnvironment(globalOperationStrategy);

		public void Deserialize(IXunitSerializationInfo info)
		{
			//globalOperationStrategy = info.GetValue<IMpiGlobalOperationStrategy>("globalOperationStrategy");
		}

		public void Serialize(IXunitSerializationInfo info)
		{
			//info.AddValue("globalOperationStrategy", globalOperationStrategy);
		}

		public override string ToString()
		{
			return "MpiEnvironment";
		}
	}
}
