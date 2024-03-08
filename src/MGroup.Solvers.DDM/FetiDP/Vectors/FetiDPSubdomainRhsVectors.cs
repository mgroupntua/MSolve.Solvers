using MGroup.LinearAlgebra.Vectors;
using MGroup.Solvers.DDM.FetiDP.Dofs;
using MGroup.Solvers.DDM.FetiDP.StiffnessMatrices;
using MGroup.Solvers.DDM.LinearSystem;

namespace MGroup.Solvers.DDM.FetiDP.Vectors
{
	public class FetiDPSubdomainRhsVectors
	{
		private readonly ISubdomainLinearSystem linearSystem;
		private readonly IFetiDPSubdomainMatrixManager matrixManagerFetiDP;
		private readonly FetiDPSubdomainDofs subdomainDofs;
		private readonly SubdomainLagranges subdomainLagranges;

		public Vector VectorFc { get; private set; }

		public Vector VectorFcCondensed { get; private set; }

		public Vector VectorFr { get; private set; }

		public Vector VectorInvKrrTimesFr { get; private set; }

		public FetiDPSubdomainRhsVectors(ISubdomainLinearSystem linearSystem, FetiDPSubdomainDofs subdomainDofs, 
			SubdomainLagranges subdomainLagranges, IFetiDPSubdomainMatrixManager matrixManagerFetiDP)
		{
			this.linearSystem = linearSystem;
			this.subdomainDofs = subdomainDofs;
			this.subdomainLagranges = subdomainLagranges;
			this.matrixManagerFetiDP = matrixManagerFetiDP;
		}

		public bool IsEmpty => VectorFr == null;

		public void CalcCondensedRhsVector()
		{
			// Static condensation: fcCondensed[s] = fc[s] - Kcr[s] * inv(Krr[s]) * fr[s]
			VectorInvKrrTimesFr = matrixManagerFetiDP.MultiplyInverseKrrTimes(VectorFr);
			VectorFcCondensed = VectorFc - matrixManagerFetiDP.MultiplyKcrTimes(VectorInvKrrTimesFr);
		}

		public void Clear()
		{
			VectorFc = null;
			VectorFcCondensed = null;
			VectorFr = null;
			VectorInvKrrTimesFr = null;
		}

		public void ExtractRhsSubvectors(Action<Vector> scaleBoundaryVector)
		{
			int[] remainderDofs = subdomainDofs.DofsRemainderToFree;
			int[] cornerDofs = subdomainDofs.DofsCornerToFree;

			Vector ff = linearSystem.RhsVector.Copy();
			scaleBoundaryVector(ff);

			this.VectorFr = ff.GetSubvector(remainderDofs);
			this.VectorFc = ff.GetSubvector(cornerDofs);
		}
	}
}
