using System;
using System.Collections.Generic;
using System.Text;
using MGroup.LinearAlgebra.Matrices;
using MGroup.Solvers.Assemblers;
using MGroup.Solvers.DDM.LinearSystem;
using MGroup.Solvers.DDM.PSM.Dofs;

namespace MGroup.Solvers.DDM.PSM.StiffnessMatrices
{
	public interface IPsmSubdomainMatrixManagerFactory<TMatrix>
		where TMatrix : class, IMatrix
	{
		ISubdomainMatrixAssembler<TMatrix> CreateAssembler();

		IPsmSubdomainMatrixManager CreateMatrixManager(
			SubdomainLinearSystem<TMatrix> linearSystem, PsmSubdomainDofs subdomainDofs);
	}
}
