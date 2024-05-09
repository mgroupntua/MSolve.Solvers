namespace MGroup.Solvers.DDM.FetiDP.CoarseProblem
{
	using MGroup.LinearAlgebra.Matrices;
	using MGroup.LinearAlgebra.Vectors;
	using MGroup.Solvers.DDM.Commons;

	public interface IFetiDPCoarseProblemGlobalMatrix
	{
		void Clear();

		void InvertGlobalScc(int numGlobalCornerDofs, IDictionary<int, int[]> subdomainToGlobalCornerDofs,
			IDictionary<int, IMatrix> subdomainMatricesScc);

		void MultiplyInverseScc(Vector input, Vector output);

		DofPermutation ReorderGlobalCornerDofs(int numGlobalCornerDofs, IDictionary<int, int[]> subdomainToGlobalCornerDofs);
	}
}
