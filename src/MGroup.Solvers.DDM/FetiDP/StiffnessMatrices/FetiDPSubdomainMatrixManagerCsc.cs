namespace MGroup.Solvers.DDM.FetiDP.StiffnessMatrices
{
	using MGroup.LinearAlgebra.Implementations;
	using MGroup.LinearAlgebra.Matrices;
	using MGroup.LinearAlgebra.SchurComplements;
	using MGroup.LinearAlgebra.SchurComplements.SubmatrixExtractors;
	using MGroup.LinearAlgebra.Triangulation;
	using MGroup.LinearAlgebra.Vectors;
	using MGroup.Solvers.Assemblers;
	using MGroup.Solvers.DDM.Commons;
	using MGroup.Solvers.DDM.FetiDP.Dofs;
	using MGroup.Solvers.DDM.LinearSystem;

	public class FetiDPSubdomainMatrixManagerCsc : IFetiDPSubdomainMatrixManager
	{
		private readonly bool clearKrrAfterFactorization;
		private readonly SubdomainLinearSystem<CsrMatrix> linearSystem;
		private readonly IImplementationProvider provider;
		private readonly FetiDPSubdomainDofs subdomainDofs;
		private readonly SubmatrixExtractorPckCsrCscSym submatrixExtractorBoundaryInternal = new SubmatrixExtractorPckCsrCscSym();
		private readonly SubmatrixExtractorFullCsrCsc submatrixExtractorCornerRemainder = new SubmatrixExtractorFullCsrCsc();

		private Matrix Kbb, Kcc;
		private CsrMatrix Kbi, Kcr;
		private CsrMatrix Kib, Krc;
		private CscMatrix Kii, Krr;
		private ILUCscFactorization inverseKii, inverseKrr;
		private DiagonalMatrix inverseKiiDiagonal;
		private Matrix Scc;

		public FetiDPSubdomainMatrixManagerCsc(IImplementationProvider provider,
			SubdomainLinearSystem<CsrMatrix> linearSystem, FetiDPSubdomainDofs subdomainDofs, bool clearKrrAfterFactorization)
		{
			this.provider = provider;
			this.linearSystem = linearSystem;
			this.subdomainDofs = subdomainDofs;
			this.clearKrrAfterFactorization = clearKrrAfterFactorization;
		}

		public bool IsEmpty => inverseKrr == null;

		public IMatrix SchurComplementOfRemainderDofs => Scc;

		public Matrix CalcInvKrrTimesKrc()
		{
			int numRemainderDofs = Krc.NumRows;
			int numCornerDofs = Krc.NumColumns;
			var result = Matrix.CreateZero(numRemainderDofs, numCornerDofs);
			for (int j = 0; j < numCornerDofs; ++j)
			{
				Vector colKrc = Krc.GetColumn(j);
				Vector colResult = inverseKrr.SolveLinearSystem(colKrc);
				result.SetSubcolumn(j, colResult, 0);
			}

			return result;
		}

		public void CalcSchurComplementOfRemainderDofs()
		{
			Scc = Matrix.CreateZero(Kcc.NumRows, Kcc.NumColumns);
			SchurComplementFullCsrCsrCsc.CalcSchurComplement(Kcc, Kcr, Krc, inverseKrr, Scc);
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
			Krc = null;
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
			Kib = null;
			Kii = null;
		}

		public void ExtractKiiKbbKib()
		{
			throw new NotImplementedException();
			//int[] boundaryRemainderToRemainder = subdomainDofs.DofsBoundaryRemainderToRemainder;
			//int[] internalToRemainder = subdomainDofs.DofsInternalToRemainder;

			//submatrixExtractorBoundaryInternal.ExtractSubmatrices(Krr, boundaryRemainderToRemainder, internalToRemainder);
			//Kbb = submatrixExtractorBoundaryInternal.Submatrix00;
			//Kbi = submatrixExtractorBoundaryInternal.Submatrix01;
			//Kib = submatrixExtractorBoundaryInternal.Submatrix10;
			//Kii = submatrixExtractorBoundaryInternal.Submatrix11;
		}

		public void ExtractKrrKccKrc()
		{
			int[] cornerToFree = subdomainDofs.DofsCornerToFree;
			int[] remainderToFree = subdomainDofs.DofsRemainderToFree;

			CsrMatrix Kff = linearSystem.Matrix;
			submatrixExtractorCornerRemainder.ExtractSubmatrices(Kff, cornerToFree, remainderToFree);
			Kcc = submatrixExtractorCornerRemainder.Submatrix00;
			Kcr = submatrixExtractorCornerRemainder.Submatrix01;
			Krc = submatrixExtractorCornerRemainder.Submatrix10;
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

				inverseKii = provider.CreateLUTriangulation();
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

			inverseKrr = provider.CreateLUTriangulation();
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

		public Vector MultiplyKibTimes(Vector vector) => Kib * vector;

		public Vector MultiplyKrcTimes(Vector vector) => Krc * vector;

		public void ReorderInternalDofs() => subdomainDofs.ReorderInternalDofs(DofPermutation.CreateNoPermutation());

		public void ReorderRemainderDofs() => subdomainDofs.ReorderRemainderDofs(DofPermutation.CreateNoPermutation());

		public class Factory : IFetiDPSubdomainMatrixManagerFactory<CsrMatrix>
		{
			private readonly bool clearKrrAfterFactorization;

			public Factory(bool clearKrrAfterFactorization = false)
			{
				this.clearKrrAfterFactorization = clearKrrAfterFactorization;
			}

			public ISubdomainMatrixAssembler<CsrMatrix> CreateAssembler() => new CsrMatrixAssembler(false);

			public IFetiDPSubdomainMatrixManager CreateMatrixManager(IImplementationProvider provider,
				SubdomainLinearSystem<CsrMatrix> linearSystem, FetiDPSubdomainDofs subdomainDofs)
			{
				return new FetiDPSubdomainMatrixManagerCsc(
					provider, linearSystem, subdomainDofs, clearKrrAfterFactorization);
			}
		}
	}
}
