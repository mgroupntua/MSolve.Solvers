namespace MGroup.Solvers.DDM.PSM.StiffnessMatrices
{
	using MGroup.LinearAlgebra.Matrices;
	using MGroup.LinearAlgebra.Reordering;
	using MGroup.LinearAlgebra.SchurComplements;
	using MGroup.LinearAlgebra.SchurComplements.SubmatrixExtractors;
	using MGroup.LinearAlgebra.Triangulation;
	using MGroup.LinearAlgebra.Vectors;
	using MGroup.Solvers.Assemblers;
	using MGroup.Solvers.DDM.Commons;
	using MGroup.Solvers.DDM.LinearSystem;
	using MGroup.Solvers.DDM.PSM.Dofs;

	public class PsmSubdomainMatrixManagerSymmetricCSparse : IPsmSubdomainMatrixManager
	{
		private readonly SubdomainLinearSystem<SymmetricCscMatrix> linearSystem;
		private readonly PsmSubdomainDofs subdomainDofs;
		private readonly OrderingAmdCSparseNet reordering = new OrderingAmdCSparseNet();
		private readonly SubmatrixExtractorCsrCscSym submatrixExtractor = new SubmatrixExtractorCsrCscSym();

		private CsrMatrix Kbb;
		private CsrMatrix Kbi;
		private SymmetricCscMatrix Kii;
		private CholeskyCSparseNet inverseKii;

		public PsmSubdomainMatrixManagerSymmetricCSparse( 
			SubdomainLinearSystem<SymmetricCscMatrix> linearSystem, PsmSubdomainDofs subdomainDofs)
		{
			this.linearSystem = linearSystem;
			this.subdomainDofs = subdomainDofs;
		}

		public bool IsEmpty => inverseKii == null;

		public IMatrixView CalcSchurComplement()
		{
			//TODO: Implement SchurComplement with A11 being in CSR format.
			TriangularUpper kbbUpper = Kbb.ExtractUpperAndDiagonalToPacked();
			var kbbSymm = SymmetricMatrix.CreateFromPackedColumnMajorArray(kbbUpper.RawData);
			return SchurComplementPckCsrSymCsc.CalcSchurComplement(kbbSymm, Kbi, inverseKii);
		}

		public void ClearSubMatrices()
		{
			Kbb = null;
			Kbi = null;
			Kii = null;
			inverseKii = null;
		}

		//TODO: Optimize this method. It is too slow.
		public void ExtractKiiKbbKib()
		{
			int[] boundaryDofs = subdomainDofs.DofsBoundaryToFree;
			int[] internalDofs = subdomainDofs.DofsInternalToFree;

			SymmetricCscMatrix Kff = linearSystem.Matrix;
			submatrixExtractor.ExtractSubmatrices(Kff, boundaryDofs, internalDofs);
			Kbb = submatrixExtractor.Submatrix00;
			Kbi = submatrixExtractor.Submatrix01;
			Kii = submatrixExtractor.Submatrix11;
		}

		public void HandleDofsWereModified()
		{
			ClearSubMatrices();
			submatrixExtractor.Clear();
		}

		public void InvertKii()
		{
			inverseKii = CholeskyCSparseNet.Factorize(Kii);
			Kii = null; // This memory is not overwritten, but it is not needed anymore either.
		}

		public Vector MultiplyInverseKii(Vector vector) => inverseKii.SolveLinearSystem(vector);

		public Vector MultiplyKbb(Vector vector) => Kbb * vector;

		public Vector MultiplyKbi(Vector vector) => Kbi.Multiply(vector, false);

		public Vector MultiplyKib(Vector vector) => Kbi.Multiply(vector, true);

		public void ReorderInternalDofs()
		{
			int[] internalDofs = subdomainDofs.DofsInternalToFree;
			SymmetricCscMatrix Kff = linearSystem.Matrix;
			(int[] rowIndicesKii, int[] colOffsetsKii) = submatrixExtractor.ExtractSparsityPattern(Kff, internalDofs);
			(int[] permutation, bool oldToNew) = reordering.FindPermutation(
				internalDofs.Length, rowIndicesKii, colOffsetsKii);

			subdomainDofs.ReorderInternalDofs(DofPermutation.Create(permutation, oldToNew));
		}

		public class Factory : IPsmSubdomainMatrixManagerFactory<SymmetricCscMatrix>
		{
			public ISubdomainMatrixAssembler<SymmetricCscMatrix> CreateAssembler() => new SymmetricCscMatrixAssembler(true);

			public IPsmSubdomainMatrixManager CreateMatrixManager(
				SubdomainLinearSystem<SymmetricCscMatrix> linearSystem, PsmSubdomainDofs subdomainDofs)
				=> new PsmSubdomainMatrixManagerSymmetricCSparse(linearSystem, subdomainDofs);
		}
	}
}
