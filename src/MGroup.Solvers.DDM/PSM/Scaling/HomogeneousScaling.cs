using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using MGroup.LinearAlgebra.Vectors;
using System.Collections.Concurrent;
using MGroup.Environments;
using MGroup.Solvers.DDM.PSM.Dofs;
using MGroup.LinearAlgebra.Distributed.Overlapping;
using System.Diagnostics;
using MGroup.LinearAlgebra.Matrices;
using MGroup.Solvers.DDM.LinearSystem;
using MGroup.Solvers.DDM.PSM.Reanalysis;

namespace MGroup.Solvers.DDM.PSM.Scaling
{
	public class HomogeneousScaling : IBoundaryDofScaling
	{
		private readonly IComputeEnvironment environment;
		private readonly Func<int, PsmSubdomainDofs> getSubdomainDofs;
		private readonly PsmReanalysisOptions reanalysis;
		private readonly ConcurrentDictionary<int, double[]> inverseMultiplicities = new ConcurrentDictionary<int, double[]>();

		public HomogeneousScaling(IComputeEnvironment environment, Func<int, PsmSubdomainDofs> getSubdomainDofs, 
			PsmReanalysisOptions reanalysis)
		{
			this.environment = environment;
			this.getSubdomainDofs = getSubdomainDofs;
			this.reanalysis = reanalysis;
		}

		public IDictionary<int, DiagonalMatrix> SubdomainMatricesWb { get; } = new ConcurrentDictionary<int, DiagonalMatrix>();

		public void CalcScalingMatrices(DistributedOverlappingIndexer boundaryDofIndexer)
		{
			bool isFirstAnalysis = inverseMultiplicities.Count == 0;
			Action<int> calcSubdomainScaling = subdomainID =>
			{
				if (isFirstAnalysis || !reanalysis.RhsVectors 
					|| reanalysis.ModifiedSubdomains.IsConnectivityModified(subdomainID))
				{
					#region log
					//Console.WriteLine($"Processing inverse multiplicities of subdomain {subdomainID}");
					//Debug.WriteLine($"Processing inverse multiplicities of subdomain {subdomainID}");
					#endregion

					int numBoundaryDofs = getSubdomainDofs(subdomainID).DofsBoundaryToFree.Length;

					var subdomainW = new double[numBoundaryDofs];
					double[] inverseMultiplicities = boundaryDofIndexer.GetLocalComponent(subdomainID).InverseMultiplicities;
					Array.Copy(inverseMultiplicities, subdomainW, numBoundaryDofs);

					this.inverseMultiplicities[subdomainID] = subdomainW;
					SubdomainMatricesWb[subdomainID] = DiagonalMatrix.CreateFromArray(subdomainW);
				}
			};
			environment.DoPerNode(calcSubdomainScaling);
		}

		public void ScaleBoundaryRhsVector(int subdomainID, Vector boundaryRhsVector)
		{
			double[] coefficients = inverseMultiplicities[subdomainID];
			Debug.Assert(boundaryRhsVector.Length == coefficients.Length);
			for (int i = 0; i < coefficients.Length; i++)
			{
				boundaryRhsVector[i] *= coefficients[i];
			}
		}
	}
}
