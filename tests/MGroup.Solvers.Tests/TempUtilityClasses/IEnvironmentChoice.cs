namespace MGroup.Solvers.Tests.TempUtilityClasses
{
	using System;
	using System.Collections.Generic;
	using System.Linq;
	using System.Text;
	using System.Threading.Tasks;

	using MGroup.Environments;
	using Xunit.Abstractions;

	public interface IEnvironmentChoice : IXunitSerializable
	{
		public IComputeEnvironment Activate();
	}
}
