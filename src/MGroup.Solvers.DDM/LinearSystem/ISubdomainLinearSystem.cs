namespace MGroup.Solvers.DDM.LinearSystem
{
	using MGroup.LinearAlgebra.Matrices;
	using MGroup.LinearAlgebra.Vectors;
	using MGroup.Solvers.DofOrdering;

	public interface ISubdomainLinearSystem
	{
		ISubdomainFreeDofOrdering DofOrdering { get; }

		IMatrix Matrix { get; }

		Vector RhsVector { get; }

		Vector Solution { get; set; }

		int SubdomainID { get; }
	}
}
