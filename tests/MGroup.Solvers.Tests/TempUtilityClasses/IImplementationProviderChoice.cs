namespace MGroup.Solvers.Tests.TempUtilityClasses
{
	using System;
	using System.Collections.Generic;
	using System.Linq;
	using System.Text;
	using System.Threading.Tasks;

	using MGroup.LinearAlgebra.Implementations;

	using Xunit.Abstractions;

	public interface IImplementationProviderChoice : IXunitSerializable
	{
		public IImplementationProvider Activate();
	}
}
