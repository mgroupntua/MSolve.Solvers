namespace MGroup.Solvers.DDM.FetiDP.CoarseProblem
{
	using MGroup.LinearAlgebra.Matrices;
	using MGroup.LinearAlgebra.Reordering;
	using MGroup.LinearAlgebra.Triangulation;
	using MGroup.LinearAlgebra.Vectors;
	using MGroup.Solvers.DDM.Commons;
	using MGroup.Solvers.DDM.SolversExtensions.Assemblers;

	public class FetiDPCoarseProblemMatrixSymmetricCSparse : IFetiDPCoarseProblemGlobalMatrix
	{
		private readonly SymmetricCscMatrixAssembler assembler = new SymmetricCscMatrixAssembler(true);
		private readonly OrderingAmdCSparseNet reordering = new OrderingAmdCSparseNet();

		private CholeskyCSparseNet inverseSccGlobal;

		public void Clear()
		{
			inverseSccGlobal = null;
			assembler.HandleDofOrderingWasModified();
		}

		public void InvertGlobalScc(int numGlobalCornerDofs, IDictionary<int, int[]> subdomainToGlobalCornerDofs, 
			IDictionary<int, IMatrix>  subdomainMatricesScc)
		{
			SymmetricCscMatrix globalScc = 
				assembler.BuildGlobalMatrix(numGlobalCornerDofs, subdomainToGlobalCornerDofs, subdomainMatricesScc);
			inverseSccGlobal = CholeskyCSparseNet.Factorize(globalScc);
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
