using System;
using System.Collections.Generic;
using System.Text;
using MGroup.LinearAlgebra.Vectors;
using MGroup.LinearAlgebra.Distributed.Overlapping;
using MGroup.LinearAlgebra.Matrices;
using MGroup.Solvers.DDM.LinearSystem;

namespace MGroup.Solvers.DDM.PSM.Scaling
{
	public interface IBoundaryDofScaling
	{
		IDictionary<int, DiagonalMatrix> SubdomainMatricesWb { get; }

		void CalcScalingMatrices(DistributedOverlappingIndexer indexer);

		void ScaleBoundaryRhsVector(int subdomainID, Vector subdomainForces); 
	}
}
