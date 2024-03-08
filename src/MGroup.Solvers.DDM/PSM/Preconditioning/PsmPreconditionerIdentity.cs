using System;
using System.Collections.Generic;
using System.Text;
using MGroup.LinearAlgebra.Matrices;
using MGroup.LinearAlgebra.Vectors;
using MGroup.Environments;
using MGroup.LinearAlgebra.Distributed.IterativeMethods.Preconditioning;
using MGroup.LinearAlgebra.Distributed.Overlapping;
using MGroup.Solvers.DDM.PSM.InterfaceProblem;

namespace MGroup.Solvers.DDM.PSM.Preconditioning
{
	public class PsmPreconditionerIdentity : IPsmPreconditioner
	{

		public PsmPreconditionerIdentity()
		{
		}

		public IPreconditioner Preconditioner { get; private set; }

		public void Calculate(IComputeEnvironment environment, DistributedOverlappingIndexer boundaryDofIndexer,
			IPsmInterfaceProblemMatrix interfaceProblemMatrix) 
		{
			Preconditioner = new IdentityPreconditioner();
		}
	}
}
