using System;
using System.Collections.Generic;
using System.Text;
using MGroup.LinearAlgebra.Iterative;
using MGroup.LinearAlgebra.Distributed.Overlapping;
using MGroup.Solvers.DDM.LinearSystem;

namespace MGroup.Solvers.DDM.PSM.InterfaceProblem
{
	public interface IPsmInterfaceProblemMatrix
	{
		DistributedOverlappingTransformation Matrix { get; }

		void Calculate(DistributedOverlappingIndexer indexer);

		double[] ExtractDiagonal(int subdomainID);
	}
}
