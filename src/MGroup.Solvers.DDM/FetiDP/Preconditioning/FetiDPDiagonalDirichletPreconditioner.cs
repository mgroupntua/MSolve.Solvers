using MGroup.Environments;
using MGroup.LinearAlgebra.Distributed.Overlapping;
using MGroup.LinearAlgebra.Matrices;
using MGroup.LinearAlgebra.Matrices.Operators;
using MGroup.LinearAlgebra.Vectors;
using MGroup.MSolve.Solution.LinearSystem;
using MGroup.Solvers.DDM.FetiDP.Dofs;
using MGroup.Solvers.DDM.FetiDP.Scaling;
using MGroup.Solvers.DDM.FetiDP.StiffnessMatrices;

namespace MGroup.Solvers.DDM.FetiDP.Preconditioning
{
	public class FetiDPDiagonalDirichletPreconditioner : IFetiDPPreconditioner
	{
		private IComputeEnvironment environment;
		private DistributedOverlappingIndexer lagrangeVectorIndexer;
		private Func<int, SubdomainLagranges> getSubdomainLagranges;
		private Func<int, IFetiDPSubdomainMatrixManager> getSubdomainMatrices;
		private IFetiDPScaling scaling;

		public void Apply(IGlobalVector input, IGlobalVector output)
		{
			DistributedOverlappingVector ye = lagrangeVectorIndexer.CheckCompatibleVector(input);
			DistributedOverlappingVector xe = lagrangeVectorIndexer.CheckCompatibleVector(output);
			//xe.Clear(); //TODO: Clear the existing local vectors instead of reallocating them

			environment.DoPerNode(s =>
			{
				IFetiDPSubdomainMatrixManager fetiDPMatrices = getSubdomainMatrices(s);
				SignedBooleanMatrixRowMajor Dbr = getSubdomainLagranges(s).MatrixDbr;
				DiagonalMatrix Wbr = scaling.SubdomainMatricesWbr[s];

				// Dbr * Wbr * (Kbb - Kbi * invKii * Kib) * Wbr * Dbr^T * y
				Vector ys = ye.LocalVectors[s];
				Vector v1 = Wbr.Multiply(Dbr.Multiply(ys, true));
				Vector v2 = fetiDPMatrices.MultiplyKbbTimes(v1);
				Vector v3 = fetiDPMatrices.MultiplyKibTimes(v1);
				v3 = fetiDPMatrices.MultiplyInverseKiiTimes(v3, true);
				v3 = fetiDPMatrices.MultiplyKbiTimes(v3);
				v2.AddIntoThis(v3);
				xe.LocalVectors[s] = Dbr.Multiply(Wbr.Multiply(v2));
			});
			xe.SumOverlappingEntries();
		}

		public void CalcSubdomainMatrices(int subdomainID)
		{
			IFetiDPSubdomainMatrixManager subdomainMatrices = getSubdomainMatrices(subdomainID);
			subdomainMatrices.ExtractKiiKbbKib();
			subdomainMatrices.InvertKii(true);
		}

		public void Initialize(IComputeEnvironment environment, DistributedOverlappingIndexer lagrangeVectorIndexer,
			Func<int, SubdomainLagranges> getSubdomainLagranges, Func<int, IFetiDPSubdomainMatrixManager> getSubdomainMatrices,
			IFetiDPScaling scaling)
		{
			this.environment = environment;
			this.lagrangeVectorIndexer = lagrangeVectorIndexer;
			this.getSubdomainLagranges = getSubdomainLagranges;
			this.getSubdomainMatrices = getSubdomainMatrices;
			this.scaling = scaling;
		}
	}
}
