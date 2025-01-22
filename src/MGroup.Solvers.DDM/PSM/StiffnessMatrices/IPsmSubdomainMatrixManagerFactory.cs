namespace MGroup.Solvers.DDM.PSM.StiffnessMatrices
{
	using MGroup.LinearAlgebra.Implementations;
	using MGroup.LinearAlgebra.Matrices;
	using MGroup.Solvers.Assemblers;
	using MGroup.Solvers.DDM.LinearSystem;
	using MGroup.Solvers.DDM.PSM.Dofs;

	public interface IPsmSubdomainMatrixManagerFactory<TMatrix>
		where TMatrix : class, IMatrix
	{
		ISubdomainMatrixAssembler<TMatrix> CreateAssembler();

		IPsmSubdomainMatrixManager CreateMatrixManager(
			IImplementationProvider provider, SubdomainLinearSystem<TMatrix> linearSystem, PsmSubdomainDofs subdomainDofs);
	}
}
