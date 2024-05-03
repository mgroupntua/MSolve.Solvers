namespace MGroup.Solvers.DDM.FetiDP.StiffnessMatrices
{
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

	public class FetiDPSubdomainMatrixManagerSymmetricCSparse : IFetiDPSubdomainMatrixManager
	{
		/// <summary>
		/// In FETI-DP Krr is also used for the preconditioner. In PFETI-DP only Krr is only used for the Schur complement of 
		/// remainder dofs
		/// </summary>
		private readonly bool clearKrrAfterFactorization;
		private readonly SubdomainLinearSystem<SymmetricCscMatrix> linearSystem;
		private readonly FetiDPSubdomainDofs subdomainDofs;
		private readonly OrderingAmdCSparseNet reordering = new OrderingAmdCSparseNet();
		private readonly SubmatrixExtractorPckCsrCscSym submatrixExtractorBoundaryInternal = new SubmatrixExtractorPckCsrCscSym();
		private readonly SubmatrixExtractorPckCsrCscSym submatrixExtractorCornerRemainder = new SubmatrixExtractorPckCsrCscSym();

		private SymmetricMatrix Kbb, Kcc;
		private CsrMatrix Kbi, Kcr;
		private SymmetricCscMatrix Kii, Krr;
		private CholeskyCSparseNet inverseKii, inverseKrr;
		private DiagonalMatrix inverseKiiDiagonal;
		private SymmetricMatrix Scc;

		public FetiDPSubdomainMatrixManagerSymmetricCSparse(SubdomainLinearSystem<SymmetricCscMatrix> linearSystem, 
			FetiDPSubdomainDofs subdomainDofs, bool clearKrrAfterFactorization)
		{
			this.linearSystem = linearSystem;
			this.subdomainDofs = subdomainDofs;
			this.clearKrrAfterFactorization = clearKrrAfterFactorization;
		}

		public bool IsEmpty => inverseKrr == null;

		public IMatrix SchurComplementOfRemainderDofs => Scc;

		public Matrix CalcInvKrrTimesKrc()
		{
			throw new NotImplementedException();
		}

		public void CalcSchurComplementOfRemainderDofs()
		{
			Scc = SymmetricMatrix.CreateZero(Kcc.Order);
			SchurComplementPckCsrSymCsc.CalcSchurComplement(Kcc, Kcr, inverseKrr, Scc);
		}

		public void ClearSubMatrices()
		{
			inverseKrr = null;
			Kcc = null;
			Kcr = null;
			Krr = null;
			Scc = null;

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
			//TODO: It would be better if these were returned by the extractor, instead of stored in its properties. 
			//		The only state that the extractor needs is its private mapping arrays
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
				inverseKii = CholeskyCSparseNet.Factorize(Kii);
			}
			Kii = null; // It has not been mutated, but it is no longer needed
		}

		public void InvertKrr()
		{
			inverseKrr = CholeskyCSparseNet.Factorize(Krr);
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
			(int[] rowIndicesKii, int[] colOffsetsKii) = submatrixExtractorBoundaryInternal.ExtractSparsityPattern(Krr, internalDofs);
			(int[] permutation, bool oldToNew) = reordering.FindPermutation(
				internalDofs.Length, rowIndicesKii, colOffsetsKii);

			subdomainDofs.ReorderInternalDofs(DofPermutation.Create(permutation, oldToNew));
		}

		public void ReorderRemainderDofs()
		{
			int[] remainderDofs = subdomainDofs.DofsRemainderToFree;
			SymmetricCscMatrix Kff = linearSystem.Matrix;
			(int[] rowIndicesKrr, int[] colOffsetsKrr) = submatrixExtractorCornerRemainder.ExtractSparsityPattern(Kff, remainderDofs);
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

			public IFetiDPSubdomainMatrixManager CreateMatrixManager(
				SubdomainLinearSystem<SymmetricCscMatrix> linearSystem, FetiDPSubdomainDofs subdomainDofs)
				=> new FetiDPSubdomainMatrixManagerSymmetricCSparse(linearSystem, subdomainDofs, clearKrrAfterFactorization);
		}
	}
}
