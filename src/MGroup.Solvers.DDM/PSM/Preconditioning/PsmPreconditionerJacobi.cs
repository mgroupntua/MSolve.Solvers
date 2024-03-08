using System;
using System.Collections.Generic;
using System.Text;
using MGroup.LinearAlgebra.Vectors;
using MGroup.Environments;
using MGroup.LinearAlgebra.Distributed;
using MGroup.LinearAlgebra.Distributed.IterativeMethods.Preconditioning;
using MGroup.LinearAlgebra.Distributed.Overlapping;
using MGroup.Solvers.DDM.PSM.InterfaceProblem;
using MGroup.Solvers.DDM.PSM.StiffnessMatrices;

namespace MGroup.Solvers.DDM.PSM.Preconditioning
{
    public class PsmPreconditionerJacobi : IPsmPreconditioner
    {
        public PsmPreconditionerJacobi()
        {
        }

        public IPreconditioner Preconditioner { get; private set; }

        public void Calculate(IComputeEnvironment environment, DistributedOverlappingIndexer boundaryDofIndexer,
            IPsmInterfaceProblemMatrix interfaceProblemMatrix)
        {
            Func<int, Vector> extractDiagonal = subdomainID 
                => Vector.CreateFromArray(interfaceProblemMatrix.ExtractDiagonal(subdomainID));
            Dictionary<int, Vector> localDiagonals = environment.CalcNodeData(extractDiagonal);
            var distributedDiagonal = new DistributedOverlappingVector(boundaryDofIndexer, localDiagonals);
            
            // All dofs belong to 2 or more subdomains and must have the stiffness contributions from all these subdomains.
            distributedDiagonal.SumOverlappingEntries();

            this.Preconditioner = new DistributedOverlappingJacobiPreconditioner(environment, distributedDiagonal);
        }
    }
}
