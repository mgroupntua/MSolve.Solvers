namespace MGroup.Solvers.DDM.FetiDP.CoarseProblem
{
	using MGroup.LinearAlgebra.Implementations;
	using MGroup.LinearAlgebra.Matrices;
	using MGroup.LinearAlgebra.Reordering;
	using MGroup.LinearAlgebra.Triangulation;
	using MGroup.LinearAlgebra.Vectors;
	using MGroup.Solvers.DDM.Commons;
	using MGroup.Solvers.DDM.SolversExtensions.Assemblers;

	public class FetiDPCoarseProblemMatrixSymmetricCsc : IFetiDPCoarseProblemGlobalMatrix
	{
		private readonly SymmetricCscMatrixAssembler assembler = new SymmetricCscMatrixAssembler(true);
		private readonly IImplementationProvider provider;
		private readonly AmdSymmetricOrdering reordering;

		private ICholeskySymmetricCsc inverseSccGlobal;

		public FetiDPCoarseProblemMatrixSymmetricCsc(IImplementationProvider provider)
		{
			this.provider = provider;
			this.reordering = new AmdSymmetricOrdering(provider);
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
			SymmetricCscMatrix globalScc =
				assembler.BuildGlobalMatrix(numGlobalCornerDofs, subdomainToGlobalCornerDofs, subdomainMatricesScc);

			if (inverseSccGlobal != null)
			{
				inverseSccGlobal.Dispose();
			}

			inverseSccGlobal = provider.CreateCholeskyTriangulation();
			inverseSccGlobal.Factorize(globalScc);
		}

		public void MultiplyInverseScc(Vector input, Vector output) => inverseSccGlobal.SolveLinearSystem(input, output);

		public DofPermutation ReorderGlobalCornerDofs(int numGlobalCornerDofs, IDictionary<int, int[]> subdomainToGlobalCornerDofs)
		{
			var pattern = SparsityPatternSymmetric.CreateEmpty(numGlobalCornerDofs);
			foreach (int s in subdomainToGlobalCornerDofs.Keys)
			{
				int[] globalDofs = subdomainToGlobalCornerDofs[s];
				pattern.ConnectIndices(globalDofs, false);
			}

			(int[] permutation, bool oldToNew) = reordering.FindPermutation(pattern);
			return DofPermutation.Create(permutation, oldToNew);
		}
	}
}
