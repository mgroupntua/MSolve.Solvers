using System;
using System.Collections.Generic;
using System.Text;
using MGroup.LinearAlgebra.Vectors;
using MGroup.MSolve.Solution.LinearSystem;
using MGroup.Solvers.DDM.LinearSystem;
using MGroup.Solvers.DDM.PSM.Dofs;
using MGroup.Solvers.DDM.PSM.StiffnessMatrices;

namespace MGroup.Solvers.DDM.PSM.Vectors
{
	public class PsmSubdomainVectors
	{
		private readonly IPsmSubdomainMatrixManager matrixManagerPsm;
		private readonly PsmSubdomainDofs subdomainDofs;
		private readonly ISubdomainLinearSystem linearSystem;

		private Vector vectorFb;
		private Vector vectorFi;

		public PsmSubdomainVectors(ISubdomainLinearSystem linearSystem, PsmSubdomainDofs subdomainDofs, 
			IPsmSubdomainMatrixManager matrixManagerPsm)
		{
			this.linearSystem = linearSystem;
			this.subdomainDofs = subdomainDofs;
			this.matrixManagerPsm = matrixManagerPsm;
		}

		public bool IsEmpty => vectorFi == null;

		public Vector CalcCondensedRhsVector()
		{
			// Static condensation: fbCondensed[s] = fb[s] - Kbi[s] * inv(Kii[s]) * fi[s]
			Vector temp = matrixManagerPsm.MultiplyInverseKii(vectorFi);
			temp = matrixManagerPsm.MultiplyKbi(temp);
			Vector fbCondensed = vectorFb - temp;

			return fbCondensed;
		}

		public void Clear()
		{
			vectorFi = null;
		}

		public void ExtractBoundaryInternalRhsVectors(Action<Vector> scaleBoundaryVector)
		{
			int[] internalDofs = subdomainDofs.DofsInternalToFree;
			int[] boundaryDofs = subdomainDofs.DofsBoundaryToFree;
			Vector ff = linearSystem.RhsVector;

			this.vectorFi = ff.GetSubvector(internalDofs);
			this.vectorFb = ff.GetSubvector(boundaryDofs);
			scaleBoundaryVector(vectorFb);
		}

		public Vector CalcSubdomainFreeSolution(Vector subdomainBoundarySolution)
		{
			int numFreeDofs = linearSystem.DofOrdering.NumFreeDofs;
			int[] boundaryDofs = subdomainDofs.DofsBoundaryToFree;
			int[] internalDofs = subdomainDofs.DofsInternalToFree;

			// ui[s] = inv(Kii[s]) * (fi[s] - Kib[s] * ub[s])
			Vector ub = subdomainBoundarySolution;
			Vector temp = matrixManagerPsm.MultiplyKib(ub);
			temp.LinearCombinationIntoThis(-1.0, vectorFi, +1);
			Vector ui = matrixManagerPsm.MultiplyInverseKii(temp);

			// Gather ub[s], ui[s] into uf[s]
			var uf = Vector.CreateZero(numFreeDofs);
			uf.CopyNonContiguouslyFrom(boundaryDofs, subdomainBoundarySolution);
			uf.CopyNonContiguouslyFrom(internalDofs, ui);

			return uf;
		}

		public void CalcStoreSubdomainFreeSolution(Vector subdomainBoundarySolution)
		{
			linearSystem.Solution = CalcSubdomainFreeSolution(subdomainBoundarySolution);
		}
	}
}
