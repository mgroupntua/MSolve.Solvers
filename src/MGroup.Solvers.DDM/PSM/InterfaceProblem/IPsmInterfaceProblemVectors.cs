using System;
using System.Collections.Generic;
using System.Text;
using MGroup.LinearAlgebra.Vectors;
using MGroup.Environments;
using MGroup.LinearAlgebra.Distributed.Overlapping;
using MGroup.Solvers.DDM.LinearSystem;
using MGroup.Solvers.DDM.PSM.Vectors;

namespace MGroup.Solvers.DDM.PSM.InterfaceProblem
{
	public interface IPsmInterfaceProblemVectors
	{
		DistributedOverlappingVector InterfaceProblemRhs { get; }

		DistributedOverlappingVector InterfaceProblemSolution { get; set; }

		// globalF = sum {Lb[s]^T * (fb[s] - Kbi[s] * inv(Kii[s]) * fi[s]) }
		void CalcInterfaceRhsVector(DistributedOverlappingIndexer indexer);

		void Clear();

	}
}
