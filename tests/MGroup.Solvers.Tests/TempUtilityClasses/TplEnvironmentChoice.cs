namespace MGroup.Solvers.Tests.TempUtilityClasses
{
	using System;
	using System.Collections.Generic;
	using System.Linq;
	using System.Text;
	using System.Threading.Tasks;

	using MGroup.Environments;
	using Xunit.Abstractions;

	public class TplEnvironmentChoice : IEnvironmentChoice
	{
		private bool optimizeBuffers;

		public TplEnvironmentChoice()
		{
			this.optimizeBuffers = false;
		}

		public TplEnvironmentChoice(bool optimizeBuffers)
		{
			this.optimizeBuffers = optimizeBuffers;
		}

		public IComputeEnvironment Activate() => new TplSharedEnvironment(optimizeBuffers);

		public void Deserialize(IXunitSerializationInfo info)
		{
			optimizeBuffers = info.GetValue<bool>("optimizeBuffers");
		}

		public void Serialize(IXunitSerializationInfo info)
		{
			info.AddValue("optimizeBuffers", optimizeBuffers);
		}

		public override string ToString()
		{
			return "TplSharedEnvironment";
		}
	}
}
