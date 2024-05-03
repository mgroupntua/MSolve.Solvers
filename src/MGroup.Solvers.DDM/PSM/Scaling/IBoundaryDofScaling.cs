namespace MGroup.Solvers.DDM.PSM.Scaling
{
	using System.Collections.Generic;

	using MGroup.LinearAlgebra.Distributed.Overlapping;
	using MGroup.LinearAlgebra.Matrices;
	using MGroup.LinearAlgebra.Vectors;

	public interface IBoundaryDofScaling
	{
		IDictionary<int, DiagonalMatrix> SubdomainMatricesWb { get; }

		void CalcScalingMatrices(DistributedOverlappingIndexer indexer);

		void ScaleBoundaryRhsVector(int subdomainID, Vector subdomainForces); 
	}
}
