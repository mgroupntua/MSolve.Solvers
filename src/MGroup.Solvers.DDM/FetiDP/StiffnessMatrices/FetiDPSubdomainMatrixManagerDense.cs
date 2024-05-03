namespace MGroup.Solvers.DDM.FetiDP.StiffnessMatrices
{
	using MGroup.LinearAlgebra.Matrices;
	using MGroup.LinearAlgebra.Vectors;
	using MGroup.Solvers.DDM.FetiDP.Dofs;
	using MGroup.Solvers.DDM.Commons;
	using MGroup.Solvers.DDM.LinearSystem;
	using MGroup.Solvers.Assemblers;

	public class FetiDPSubdomainMatrixManagerDense : IFetiDPSubdomainMatrixManager
	{
		private readonly SubdomainLinearSystem<Matrix> linearSystem;
		private readonly FetiDPSubdomainDofs subdomainDofs;

		private Matrix Kbb, Kbi, Kib, Kii;
		private Matrix Kcc, Kcr, Krc, Krr;
		private Matrix inverseKii, inverseKrr;
		private DiagonalMatrix inverseKiiDiagonal;
		private Matrix Scc;

		public FetiDPSubdomainMatrixManagerDense(SubdomainLinearSystem<Matrix> linearSystem, FetiDPSubdomainDofs subdomainDofs)
		{
			this.linearSystem = linearSystem;
			this.subdomainDofs = subdomainDofs;
		}

		public bool IsEmpty => inverseKrr == null;

		public IMatrix SchurComplementOfRemainderDofs => Scc;

		public Matrix CalcInvKrrTimesKrc()
		{
			throw new NotImplementedException();
		}

		public void CalcSchurComplementOfRemainderDofs()
		{
			Scc = Kcc - Kcr * inverseKrr * Krc;
		}

		public void ClearSubMatrices()
		{
			inverseKrr = null;
			Kcc = null;
			Kcr = null;
			Krc = null;
			Krr = null;
			Scc = null;

			inverseKii = null;
			inverseKiiDiagonal = null;
			Kbb = null;
			Kbi = null;
			Kib = null;
			Kii = null;
		}

		public void ExtractKiiKbbKib()
		{
			(int[] internalToFree, int[] boundaryRemainderToFree) = subdomainDofs.RemapDofs();
			Matrix Kff = linearSystem.Matrix;
			Kbb = Kff.GetSubmatrix(boundaryRemainderToFree, boundaryRemainderToFree);
			Kbi = Kff.GetSubmatrix(boundaryRemainderToFree, internalToFree);
			Kib = Kff.GetSubmatrix(internalToFree, boundaryRemainderToFree);
			Kii = Kff.GetSubmatrix(internalToFree, internalToFree);
		}

		public void ExtractKrrKccKrc()
		{
			int[] cornerToFree = subdomainDofs.DofsCornerToFree;
			int[] remainderToFree = subdomainDofs.DofsRemainderToFree;
			Matrix Kff = linearSystem.Matrix;
			Kcc = Kff.GetSubmatrix(cornerToFree, cornerToFree);
			Kcr = Kff.GetSubmatrix(cornerToFree, remainderToFree);
			Krc = Kff.GetSubmatrix(remainderToFree, cornerToFree);
			Krr = Kff.GetSubmatrix(remainderToFree, remainderToFree);
		}

		public void HandleDofsWereModified() => ClearSubMatrices();

		public void InvertKii(bool diagonalOnly)
		{
			if (diagonalOnly)
			{
				inverseKiiDiagonal = DiagonalMatrix.CreateFromArray(Kii.GetDiagonalAsArray());
				inverseKiiDiagonal.Invert();
			}
			else
			{
				inverseKii = Kii.Invert();
			}
			Kii = null; // Kii has been overwritten or is not needed anymore
		}

		public void InvertKrr()
		{
			inverseKrr = Krr.Invert();
			Krr = null; // Krr has been overwritten
		}

		public Vector MultiplyInverseKiiTimes(Vector vector, bool diagonalOnly)
		{
			if (diagonalOnly)
			{
				return inverseKiiDiagonal.Multiply(vector);
			}
			else
			{
				return inverseKii* vector;
			}
		}

		public Vector MultiplyInverseKrrTimes(Vector vector) => inverseKrr * vector;

		public Vector MultiplyKbbTimes(Vector vector) => Kbb * vector;

		public Vector MultiplyKbiTimes(Vector vector) => Kbi * vector;

		public Vector MultiplyKccTimes(Vector vector) => Kcc * vector;

		public Vector MultiplyKcrTimes(Vector vector) => Kcr * vector;

		public Vector MultiplyKibTimes(Vector vector) => Kib * vector;

		public Vector MultiplyKrcTimes(Vector vector) => Krc * vector;

		public void ReorderInternalDofs() => subdomainDofs.ReorderInternalDofs(DofPermutation.CreateNoPermutation());

		public void ReorderRemainderDofs() => subdomainDofs.ReorderRemainderDofs(DofPermutation.CreateNoPermutation());

		public class Factory : IFetiDPSubdomainMatrixManagerFactory<Matrix>
		{
			public ISubdomainMatrixAssembler<Matrix> CreateAssembler() => new DenseMatrixAssembler();

			public IFetiDPSubdomainMatrixManager CreateMatrixManager(
				SubdomainLinearSystem<Matrix> linearSystem, FetiDPSubdomainDofs subdomainDofs)
				=> new FetiDPSubdomainMatrixManagerDense(linearSystem, subdomainDofs);
		}
	}
}
