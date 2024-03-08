using System;
using System.Collections.Generic;
using System.Text;
using MGroup.Environments;
using MGroup.LinearAlgebra.Distributed.IterativeMethods.Preconditioning;
using MGroup.LinearAlgebra.Distributed.Overlapping;
using MGroup.Solvers.DDM.PSM.InterfaceProblem;

namespace MGroup.Solvers.DDM.PSM.Preconditioning
{
	public interface IPsmPreconditioner
	{
		void Calculate(IComputeEnvironment environment, DistributedOverlappingIndexer boundaryDofIndexer, 
			IPsmInterfaceProblemMatrix interfaceProblemMatrix);

		IPreconditioner Preconditioner { get; }
	}
}
