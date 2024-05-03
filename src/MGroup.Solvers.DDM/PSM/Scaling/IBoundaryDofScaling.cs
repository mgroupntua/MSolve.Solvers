namespace MGroup.Solvers.DDM.PSM.Scaling
{
	using System.Collections.Generic;
	using MGroup.LinearAlgebra.Vectors;
	using MGroup.LinearAlgebra.Distributed.Overlapping;
	using MGroup.LinearAlgebra.Matrices;

	public interface IBoundaryDofScaling
	{
		IDictionary<int, DiagonalMatrix> SubdomainMatricesWb { get; }

		void CalcScalingMatrices(DistributedOverlappingIndexer indexer);

		void ScaleBoundaryRhsVector(int subdomainID, Vector subdomainForces); 
	}
}
