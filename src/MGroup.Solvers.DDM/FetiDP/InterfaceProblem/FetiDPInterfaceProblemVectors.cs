using MGroup.Environments;
using MGroup.LinearAlgebra.Distributed.Overlapping;
using MGroup.LinearAlgebra.Vectors;
using MGroup.Solvers.DDM.FetiDP.CoarseProblem;
using MGroup.Solvers.DDM.FetiDP.Dofs;
using MGroup.Solvers.DDM.FetiDP.StiffnessMatrices;
using MGroup.Solvers.DDM.FetiDP.Vectors;

namespace MGroup.Solvers.DDM.FetiDP.InterfaceProblem
{
	public class FetiDPInterfaceProblemVectors : IFetiDPInterfaceProblemVectors
	{
		private readonly IFetiDPCoarseProblem coarseProblem;
		private readonly IComputeEnvironment environment;
		private readonly IDictionary<int, SubdomainLagranges> subdomainLagranges;
		private readonly IDictionary<int, IFetiDPSubdomainMatrixManager> subdomainMatrices;
		private readonly IDictionary<int, FetiDPSubdomainRhsVectors> subdomainVectors;

		public FetiDPInterfaceProblemVectors(IComputeEnvironment environment, IFetiDPCoarseProblem coarseProblem, 
			IDictionary<int, SubdomainLagranges> subdomainLagranges,
			IDictionary<int, IFetiDPSubdomainMatrixManager> subdomainMatrices, 
			IDictionary<int, FetiDPSubdomainRhsVectors> subdomainVectors)
		{
			this.environment = environment;
			this.coarseProblem = coarseProblem;
			this.subdomainLagranges = subdomainLagranges;
			this.subdomainMatrices = subdomainMatrices;
			this.subdomainVectors = subdomainVectors;
		}

		public DistributedOverlappingVector InterfaceProblemRhs { get; private set; }

		public DistributedOverlappingVector InterfaceProblemSolution { get; set; }

		public void CalcInterfaceRhsVector(DistributedOverlappingIndexer lagrangeVectorIndexer)
		{
			Dictionary<int, Vector> yce = environment.CalcNodeData(s => subdomainVectors[s].VectorFcCondensed);
			Dictionary<int, Vector> xce = environment.CalcNodeData(s => Vector.CreateZero(yce[s].Length));
			coarseProblem.SolveCoarseProblem(yce, xce);

			Dictionary<int, Vector> rhsVectors = environment.CalcNodeData(s =>
			{
				Vector v1s = subdomainVectors[s].VectorInvKrrTimesFr;
				Vector v2s = subdomainMatrices[s].MultiplyInverseKrrTimes(subdomainMatrices[s].MultiplyKrcTimes(xce[s]));
				v2s.LinearCombinationIntoThis(-1, v1s, +1);
				return subdomainLagranges[s].MatrixDr.Multiply(v2s);
			});
			InterfaceProblemRhs = new DistributedOverlappingVector(lagrangeVectorIndexer, rhsVectors);
			InterfaceProblemRhs.SumOverlappingEntries();
		}

		public void Clear()
		{
			InterfaceProblemRhs = null;
		}
	}
}
