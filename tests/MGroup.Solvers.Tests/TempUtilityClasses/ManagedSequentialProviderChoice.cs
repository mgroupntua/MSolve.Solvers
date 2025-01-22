namespace MGroup.Solvers.Tests.TempUtilityClasses
{
	using System;
	using System.Collections.Generic;
	using System.Linq;
	using System.Text;
	using System.Threading.Tasks;

	using MGroup.LinearAlgebra.Implementations;
	using MGroup.LinearAlgebra.Implementations.Managed;
	using MGroup.LinearAlgebra.Triangulation;
	using Xunit.Abstractions;

	public class ManagedSequentialProviderChoice : IImplementationProviderChoice
	{
		private double luPivotTolerance;

		public ManagedSequentialProviderChoice()
		{
			this.luPivotTolerance = LUCSparseNet.DefaultPivotTolerance;
		}

		public ManagedSequentialProviderChoice(double luPivotTolerance)
		{
			this.luPivotTolerance = luPivotTolerance;
		}

		public IImplementationProvider Activate() => new ManagedSequentialImplementationProvider(luPivotTolerance);

		public void Deserialize(IXunitSerializationInfo info)
		{
			luPivotTolerance = info.GetValue<double>("luPivotTolerance");
		}

		public void Serialize(IXunitSerializationInfo info)
		{
			info.AddValue("luPivotTolerance", luPivotTolerance);
		}

		public override string ToString()
		{
			return "ManagedSequentialImplementationProvider";
		}
	}
}
