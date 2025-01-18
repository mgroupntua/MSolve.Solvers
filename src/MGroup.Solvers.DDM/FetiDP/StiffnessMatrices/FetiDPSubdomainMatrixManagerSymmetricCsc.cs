namespace MGroup.Solvers.DDM.FetiDP.StiffnessMatrices
{
	using MGroup.LinearAlgebra.Implementations;
	using MGroup.LinearAlgebra.Matrices;
	using MGroup.LinearAlgebra.Reordering;
	using MGroup.LinearAlgebra.SchurComplements;
	using MGroup.LinearAlgebra.SchurComplements.SubmatrixExtractors;
	using MGroup.LinearAlgebra.Triangulation;
	using MGroup.LinearAlgebra.Vectors;
	using MGroup.Solvers.Assemblers;
	using MGroup.Solvers.DDM.Commons;
	using MGroup.Solvers.DDM.FetiDP.Dofs;
	using MGroup.Solvers.DDM.LinearSystem;

	public class FetiDPSubdomainMatrixManagerSymmetricCsc : IFetiDPSubdomainMatrixManager
	{
		/// <summary>
		/// In FETI-DP Krr is also used for the preconditioner. In PFETI-DP only Krr is only used for the Schur complement of 
		/// remainder dofs.
		/// </summary>
		private readonly bool clearKrrAfterFactorization;
		private readonly SubdomainLinearSystem<SymmetricCscMatrix> linearSystem;
		private readonly IImplementationProvider provider;
		private readonly AmdSymmetricOrdering reordering;
		private readonly FetiDPSubdomainDofs subdomainDofs;
		private readonly SubmatrixExtractorPckCsrCscSym submatrixExtractorBoundaryInternal = new SubmatrixExtractorPckCsrCscSym();
		private readonly SubmatrixExtractorPckCsrCscSym submatrixExtractorCornerRemainder = new SubmatrixExtractorPckCsrCscSym();

		private SymmetricMatrix Kbb, Kcc;
		private CsrMatrix Kbi, Kcr;
		private SymmetricCscMatrix Kii, Krr;
		private ICholeskySymmetricCsc inverseKii, inverseKrr;
		private DiagonalMatrix inverseKiiDiagonal;
		private SymmetricMatrix Scc;

		public FetiDPSubdomainMatrixManagerSymmetricCsc(
			IImplementationProvider provider, SubdomainLinearSystem<SymmetricCscMatrix> linearSystem, 
			FetiDPSubdomainDofs subdomainDofs, bool clearKrrAfterFactorization)
		{
			this.provider = provider;
			this.linearSystem = linearSystem;
			this.subdomainDofs = subdomainDofs;
			this.clearKrrAfterFactorization = clearKrrAfterFactorization;
			this.reordering = new AmdSymmetricOrdering(provider);
		}

		public bool IsEmpty => inverseKrr == null;

		public IMatrix SchurComplementOfRemainderDofs => Scc;

		public Matrix CalcInvKrrTimesKrc()
		{
			int numCornerDofs = Kcr.NumRows;
			int numRemainderDofs = Kcr.NumColumns;
			var result = Matrix.CreateZero(numRemainderDofs, numCornerDofs);
			for (int j = 0; j < numCornerDofs; ++j)
			{
				Vector colKrc = Kcr.GetRow(j);
				Vector colResult = inverseKrr.SolveLinearSystem(colKrc);
				result.SetSubcolumn(j, colResult, 0);
			}

			return result;
		}

		public void CalcSchurComplementOfRemainderDofs()
		{
			Scc = SymmetricMatrix.CreateZero(Kcc.Order);
			SchurComplementPckCsrSymCsc.CalcSchurComplement(Kcc, Kcr, inverseKrr, Scc);
		}

		public void ClearSubMatrices()
		{
			if (inverseKrr != null)
			{
				inverseKrr.Dispose();
			}

			inverseKrr = null;
			Kcc = null;
			Kcr = null;
			Krr = null;
			Scc = null;

			if (inverseKii != null)
			{
				inverseKii.Dispose();
			}

			inverseKii = null;
			inverseKiiDiagonal = null;
			Kbb = null;
			Kbi = null;
			Kii = null;
		}

		public void ExtractKiiKbbKib()
		{
			int[] boundaryRemainderToRemainder = subdomainDofs.DofsBoundaryRemainderToRemainder;
			int[] internalToRemainder = subdomainDofs.DofsInternalToRemainder;

			submatrixExtractorBoundaryInternal.ExtractSubmatrices(Krr, boundaryRemainderToRemainder, internalToRemainder);
			Kbb = submatrixExtractorBoundaryInternal.Submatrix00;
			Kbi = submatrixExtractorBoundaryInternal.Submatrix01;
			Kii = submatrixExtractorBoundaryInternal.Submatrix11;
		}

		public void ExtractKrrKccKrc()
		{
			int[] cornerToFree = subdomainDofs.DofsCornerToFree;
			int[] remainderToFree = subdomainDofs.DofsRemainderToFree;

			SymmetricCscMatrix Kff = linearSystem.Matrix;
			submatrixExtractorCornerRemainder.ExtractSubmatrices(Kff, cornerToFree, remainderToFree);
			Kcc = submatrixExtractorCornerRemainder.Submatrix00;
			Kcr = submatrixExtractorCornerRemainder.Submatrix01;
			Krr = submatrixExtractorCornerRemainder.Submatrix11;

			submatrixExtractorBoundaryInternal.Clear();
		}

		public void HandleDofsWereModified()
		{
			ClearSubMatrices();
			submatrixExtractorCornerRemainder.Clear();
		}

		public void InvertKii(bool diagonalOnly)
		{
			if (diagonalOnly)
			{
				inverseKiiDiagonal = DiagonalMatrix.CreateFromArray(((IDiagonalAccessible)Kii).GetDiagonalAsArray());
				inverseKiiDiagonal.Invert();
			}
			else
			{
				if (inverseKii != null)
				{
					inverseKii.Dispose();
				}

				inverseKii = provider.CreateSymmetricCscTriangulation(true);
				inverseKii.Factorize(Kii);
			}

			Kii = null; // It has not been mutated, but it is no longer needed
		}

		public void InvertKrr()
		{
			if (inverseKrr != null)
			{
				inverseKrr.Dispose();
			}

			inverseKrr = provider.CreateSymmetricCscTriangulation(true);
			inverseKrr.Factorize(Krr);

			if (clearKrrAfterFactorization)
			{
				Krr = null; // It has not been mutated, but it is no longer needed
			}
		}

		public Vector MultiplyInverseKiiTimes(Vector vector, bool diagonalOnly)
		{
			if (diagonalOnly)
			{
				return inverseKiiDiagonal.Multiply(vector);
			}
			else
			{
				return inverseKii.SolveLinearSystem(vector);
			}
		}

		public Vector MultiplyInverseKrrTimes(Vector vector) => inverseKrr.SolveLinearSystem(vector);

		public Vector MultiplyKbbTimes(Vector vector) => Kbb * vector;

		public Vector MultiplyKbiTimes(Vector vector) => Kbi * vector;

		public Vector MultiplyKccTimes(Vector vector) => Kcc * vector;

		public Vector MultiplyKcrTimes(Vector vector) => Kcr * vector;

		public Vector MultiplyKibTimes(Vector vector) => Kbi.Multiply(vector, true);

		public Vector MultiplyKrcTimes(Vector vector) => Kcr.Multiply(vector, true);

		public void ReorderInternalDofs()
		{
			int[] internalDofs = subdomainDofs.DofsInternalToRemainder;
			(int[] rowIndicesKii, int[] colOffsetsKii) =
				submatrixExtractorBoundaryInternal.ExtractSparsityPattern(Krr, internalDofs);
			(int[] permutation, bool oldToNew) = reordering.FindPermutation(
				internalDofs.Length, rowIndicesKii, colOffsetsKii);

			subdomainDofs.ReorderInternalDofs(DofPermutation.Create(permutation, oldToNew));
		}

		public void ReorderRemainderDofs()
		{
			int[] remainderDofs = subdomainDofs.DofsRemainderToFree;
			SymmetricCscMatrix Kff = linearSystem.Matrix;
			(int[] rowIndicesKrr, int[] colOffsetsKrr) =
				submatrixExtractorCornerRemainder.ExtractSparsityPattern(Kff, remainderDofs);
			(int[] permutation, bool oldToNew) = reordering.FindPermutation(
				remainderDofs.Length, rowIndicesKrr, colOffsetsKrr);

			subdomainDofs.ReorderRemainderDofs(DofPermutation.Create(permutation, oldToNew));
		}

		public class Factory : IFetiDPSubdomainMatrixManagerFactory<SymmetricCscMatrix>
		{
			private readonly bool clearKrrAfterFactorization;

			public Factory(bool clearKrrAfterFactorization = false)
			{
				this.clearKrrAfterFactorization = clearKrrAfterFactorization;
			}

			public ISubdomainMatrixAssembler<SymmetricCscMatrix> CreateAssembler() => new SymmetricCscMatrixAssembler(true);

			public IFetiDPSubdomainMatrixManager CreateMatrixManager(IImplementationProvider provider,
				SubdomainLinearSystem<SymmetricCscMatrix> linearSystem, FetiDPSubdomainDofs subdomainDofs)
			{
				return new FetiDPSubdomainMatrixManagerSymmetricCsc(
					provider, linearSystem, subdomainDofs, clearKrrAfterFactorization);
			}
		}
	}
}
