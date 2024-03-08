using MGroup.LinearAlgebra.Matrices;
using MGroup.Solvers.Assemblers;
using MGroup.Solvers.DDM.FetiDP.Dofs;
using MGroup.Solvers.DDM.LinearSystem;

namespace MGroup.Solvers.DDM.FetiDP.StiffnessMatrices
{
	public interface IFetiDPSubdomainMatrixManagerFactory<TMatrix>
		where TMatrix : class, IMatrix
	{
		ISubdomainMatrixAssembler<TMatrix> CreateAssembler();

		IFetiDPSubdomainMatrixManager CreateMatrixManager(
			SubdomainLinearSystem<TMatrix> linearSystem, FetiDPSubdomainDofs subdomainDofs);
	}
}
