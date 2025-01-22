namespace MGroup.Solvers.Tests.TempUtilityClasses
{
	using System;
	using System.Collections.Generic;
	using System.Linq;
	using System.Text;
	using System.Threading.Tasks;

	using MGroup.LinearAlgebra.Implementations;
	using MGroup.LinearAlgebra.Implementations.NativeWin64;
	using MGroup.LinearAlgebra.Triangulation;
	using Xunit.Abstractions;

	public class NativeWin64ProviderChoice : IImplementationProviderChoice
	{
		private double luPivotTolerance;
		private bool superNodalCholesky;

		public NativeWin64ProviderChoice()
		{
			this.luPivotTolerance = LUCSparseNet.DefaultPivotTolerance;
			this.superNodalCholesky = true;
		}

		public NativeWin64ProviderChoice(double luPivotTolerance, bool superNodalCholesky)
		{
			this.luPivotTolerance = luPivotTolerance;
			this.superNodalCholesky = superNodalCholesky;
		}

		public IImplementationProvider Activate()
			=> new NativeWin64ImplementationProvider(luPivotTolerance, superNodalCholesky);

		public void Deserialize(IXunitSerializationInfo info)
		{
			luPivotTolerance = info.GetValue<double>("luPivotTolerance");
			superNodalCholesky = info.GetValue<bool>("superNodalCholesky");
		}

		public void Serialize(IXunitSerializationInfo info)
		{
			info.AddValue("luPivotTolerance", luPivotTolerance);
			info.AddValue("superNodalCholesky", superNodalCholesky);
		}

		public override string ToString()
		{
			return "NativeWin64ImplementationProvider";
		}
	}
}
