namespace MGroup.Solvers.DDM.FetiDP.Scaling
{
	using MGroup.LinearAlgebra.Matrices;
	using MGroup.LinearAlgebra.Vectors;

	public interface IFetiDPScaling
	{
		IDictionary<int, DiagonalMatrix> SubdomainMatricesWbr { get; }

		void CalcScalingMatrices();

		void ScaleSubdomainRhsVector(int subdomainID, Vector rhsAtFreeDofs);

		void ScaleSubdomainSolutionVector(int subdomainID, Vector solutionAtFreeDofs);
	}
}
