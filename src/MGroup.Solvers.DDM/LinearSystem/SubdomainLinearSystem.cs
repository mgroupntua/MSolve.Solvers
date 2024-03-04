using System;
using System.Collections.Generic;
using System.Text;
using MGroup.LinearAlgebra.Matrices;
using MGroup.LinearAlgebra.Vectors;
using MGroup.Solvers.DofOrdering;

namespace MGroup.Solvers.DDM.LinearSystem
{
	public class SubdomainLinearSystem<TMatrix> : ISubdomainLinearSystem
		where TMatrix : class, IMatrix
	{
		private readonly DistributedAlgebraicModel<TMatrix> algebraicModel;

		public SubdomainLinearSystem(DistributedAlgebraicModel<TMatrix> algebraicModel, int subdomainID)
		{
			this.SubdomainID = subdomainID;
			this.algebraicModel = algebraicModel;
		}

		public ISubdomainFreeDofOrdering DofOrdering => algebraicModel.SubdomainFreeDofOrderings[SubdomainID];

		IMatrix ISubdomainLinearSystem.Matrix => this.Matrix;

		public TMatrix Matrix
		{
			get => algebraicModel.LinearSystem.Matrix.LocalMatrices[SubdomainID];
		}

		public Vector RhsVector 
		{
			get => algebraicModel.LinearSystem.RhsVector.LocalVectors[SubdomainID];
		}

		public Vector Solution 
		{
			get => algebraicModel.LinearSystem.Solution.LocalVectors[SubdomainID];
			set => algebraicModel.LinearSystem.Solution.LocalVectors[SubdomainID] = value;
		}

		public int SubdomainID { get; }
	}
}
