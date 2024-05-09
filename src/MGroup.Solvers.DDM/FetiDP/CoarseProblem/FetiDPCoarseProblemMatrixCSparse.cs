namespace MGroup.Solvers.DDM.FetiDP.CoarseProblem
{
	using MGroup.LinearAlgebra.Matrices;
	using MGroup.LinearAlgebra.Triangulation;
	using MGroup.LinearAlgebra.Vectors;
	using MGroup.Solvers.DDM.Commons;
	using MGroup.Solvers.DDM.SolversExtensions.Assemblers;

	public class FetiDPCoarseProblemMatrixCSparse : IFetiDPCoarseProblemGlobalMatrix
	{
		private readonly CscMatrixAssembler assembler = new CscMatrixAssembler(false, true);
		private LUCSparseNet inverseSccGlobal;

		public void Clear() 
		{
			inverseSccGlobal = null;
			assembler.HandleDofOrderingWasModified();
		}

		public void InvertGlobalScc(int numGlobalCornerDofs, IDictionary<int, int[]> subdomainToGlobalCornerDofs, 
			IDictionary<int, IMatrix>  subdomainMatricesScc)
		{
			CscMatrix globalScc = 
				assembler.BuildGlobalMatrix(numGlobalCornerDofs, subdomainToGlobalCornerDofs, subdomainMatricesScc);
			inverseSccGlobal = LUCSparseNet.Factorize(globalScc);
		}

		public void MultiplyInverseScc(Vector input, Vector output) => inverseSccGlobal.SolveLinearSystem(input, output);

		public DofPermutation ReorderGlobalCornerDofs(int numGlobalCornerDofs, IDictionary<int, int[]> subdomainToGlobalCornerDofs)
			=> DofPermutation.CreateNoPermutation();
	}
}
