namespace MGroup.Solvers.DDM.FetiDP.CoarseProblem
{
	using MGroup.LinearAlgebra.Implementations;
	using MGroup.LinearAlgebra.Matrices;
	using MGroup.LinearAlgebra.Triangulation;
	using MGroup.LinearAlgebra.Vectors;
	using MGroup.Solvers.DDM.Commons;
	using MGroup.Solvers.DDM.SolversExtensions.Assemblers;

	public class FetiDPCoarseProblemMatrixCsc : IFetiDPCoarseProblemGlobalMatrix
	{
		private readonly CscMatrixAssembler assembler = new CscMatrixAssembler(false, true);
		private readonly IImplementationProvider provider;

		private ILUCscFactorization inverseSccGlobal;

		public FetiDPCoarseProblemMatrixCsc(IImplementationProvider provider)
		{
			this.provider = provider;
		}

		public void Clear()
		{
			if (inverseSccGlobal != null)
			{
				inverseSccGlobal.Dispose();
			}

			inverseSccGlobal = null;
			assembler.HandleDofOrderingWasModified();
		}

		public void InvertGlobalScc(int numGlobalCornerDofs, IDictionary<int, int[]> subdomainToGlobalCornerDofs,
			IDictionary<int, IMatrix>  subdomainMatricesScc)
		{
			CscMatrix globalScc =
				assembler.BuildGlobalMatrix(numGlobalCornerDofs, subdomainToGlobalCornerDofs, subdomainMatricesScc);

			if (inverseSccGlobal != null)
			{
				inverseSccGlobal.Dispose();
			}

			inverseSccGlobal = provider.CreateLUTriangulation();
			inverseSccGlobal.Factorize(globalScc);
		}

		public void MultiplyInverseScc(Vector input, Vector output) => inverseSccGlobal.SolveLinearSystem(input, output);

		public DofPermutation ReorderGlobalCornerDofs(
			int numGlobalCornerDofs, IDictionary<int, int[]> subdomainToGlobalCornerDofs)
			=> DofPermutation.CreateNoPermutation();
	}
}
