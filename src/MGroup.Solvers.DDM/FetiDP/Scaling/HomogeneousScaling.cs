using System.Collections.Concurrent;

using MGroup.Environments;
using MGroup.LinearAlgebra.Matrices;
using MGroup.LinearAlgebra.Vectors;
using MGroup.MSolve.Discretization.Entities;
using MGroup.Solvers.DDM.FetiDP.Dofs;
using MGroup.Solvers.DDM.LagrangeMultipliers;

namespace MGroup.Solvers.DDM.FetiDP.Scaling
{
	public class HomogeneousScaling : IFetiDPScaling
	{
		private readonly IComputeEnvironment environment;
		private readonly IModel model;
		private readonly Func<int, FetiDPSubdomainDofs> getSubdomainDofs;
		private readonly FetiDPReanalysisOptions reanalysis;
		private readonly ConcurrentDictionary<int, double[]> inverseMultiplicitiesBoundaryRemainder 
			= new ConcurrentDictionary<int, double[]>();
		private readonly ConcurrentDictionary<int, double[]> inverseMultiplicitiesCorner
			= new ConcurrentDictionary<int, double[]>();

		public HomogeneousScaling(IComputeEnvironment environment, IModel model, 
			Func<int, FetiDPSubdomainDofs> getSubdomainDofs, ICrossPointStrategy crossPointStrategy, 
			FetiDPReanalysisOptions reanalysis)
		{
			this.environment = environment;
			this.model = model;
			this.getSubdomainDofs = getSubdomainDofs;
			this.reanalysis = reanalysis;

			if (!(crossPointStrategy is FullyRedundantLagranges))
			{
				throw new NotImplementedException();
			}
		}

		public IDictionary<int, DiagonalMatrix> SubdomainMatricesWbr { get; } = new ConcurrentDictionary<int, DiagonalMatrix>();

		public void CalcScalingMatrices()
		{
			bool isFirstAnalysis = inverseMultiplicitiesBoundaryRemainder.Count == 0;
			Action<int> calcSubdomainScaling = subdomainID =>
			{
				if (isFirstAnalysis || !reanalysis.RhsVectors
					|| reanalysis.ModifiedSubdomains.IsConnectivityModified(subdomainID))
				{
					#region log
					//Console.WriteLine($"Processing inverse multiplicities of subdomain {subdomainID}");
					//Debug.WriteLine($"Processing inverse multiplicities of subdomain {subdomainID}");
					#endregion

					FetiDPSubdomainDofs dofs = getSubdomainDofs(subdomainID);
					var Wbr = new double[dofs.DofsBoundaryRemainderToRemainder.Length];
					foreach ((int nodeID, int dofID, int brIndex) in dofs.DofOrderingBoundaryRemainder)
					{
						int multiplicity = model.GetNode(nodeID).Subdomains.Count;
						Wbr[brIndex] = 1.0 / multiplicity;
					}

					var Wc = new double[dofs.DofsCornerToFree.Length];
					{
						foreach ((int nodeID, int dofID, int cIndex) in dofs.DofOrderingCorner)
						{
							int multiplicity = model.GetNode(nodeID).Subdomains.Count;
							Wc[cIndex] = 1.0 / multiplicity;
						}
					}

					this.SubdomainMatricesWbr[subdomainID] = DiagonalMatrix.CreateFromArray(Wbr);
					this.inverseMultiplicitiesBoundaryRemainder[subdomainID] = Wbr;
					this.inverseMultiplicitiesCorner[subdomainID] = Wc;
				}
			};
			environment.DoPerNode(calcSubdomainScaling);
		}

		public void ScaleSubdomainRhsVector(int subdomainID, Vector rhsAtFreeDofs)
		{
			FetiDPSubdomainDofs subdomainDofs = getSubdomainDofs(subdomainID);
			int[] cornerToFree = subdomainDofs.DofsCornerToFree;
			int[] remainderToFree = subdomainDofs.DofsRemainderToFree;
			int[] boundaryRemainderToRemainder = subdomainDofs.DofsBoundaryRemainderToRemainder;

			double[] Wc = inverseMultiplicitiesCorner[subdomainID];
			for (int cornerDofIdx = 0; cornerDofIdx < Wc.Length; ++cornerDofIdx)
			{
				int freeDofIdx = cornerToFree[cornerDofIdx];
				rhsAtFreeDofs[freeDofIdx] *= Wc[cornerDofIdx];
			}

			double[] Wbr = inverseMultiplicitiesBoundaryRemainder[subdomainID];
			for (int brDofIdx = 0; brDofIdx < Wbr.Length; ++brDofIdx)
			{
				int freeDofIdx = remainderToFree[boundaryRemainderToRemainder[brDofIdx]];
				rhsAtFreeDofs[freeDofIdx] *= Wbr[brDofIdx];
			}
		}

		public void ScaleSubdomainSolutionVector(int subdomainID, Vector solutionAtFreeDofs) 
			=> ScaleSubdomainRhsVector(subdomainID, solutionAtFreeDofs);
	}
}
