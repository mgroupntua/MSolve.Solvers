using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using MGroup.LinearAlgebra.Matrices;
using MGroup.LinearAlgebra.Vectors;
using System.Collections.Concurrent;
using MGroup.Environments;
using MGroup.Solvers.DDM.PSM.Dofs;
using MGroup.LinearAlgebra.Distributed.Overlapping;
using MGroup.MSolve.Discretization;
using MGroup.MSolve.DataStructures;
using MGroup.Solvers.DDM.LinearSystem;
using System.Diagnostics;

//TODO: If Jacobi preconditioning is used, then the most time consuming part of finding the relative stiffnesses 
//		(distributedVector.SumOverlappingEntries()) is done there too. The two objects should synchronize to only do that once. 
namespace MGroup.Solvers.DDM.PSM.Scaling
{
	public class HeterogeneousScaling : IBoundaryDofScaling
	{
		private readonly IComputeEnvironment environment;
		private readonly ISubdomainTopology subdomainTopology;
		private readonly Func<int, ISubdomainLinearSystem> getSubdomainLinearSystem;
		private readonly Func<int, PsmSubdomainDofs> getSubdomainDofs;

		private readonly ConcurrentDictionary<int, double[]> relativeStiffnesses = new ConcurrentDictionary<int, double[]>();

		public HeterogeneousScaling(IComputeEnvironment environment, ISubdomainTopology subdomainTopology,
			Func<int, ISubdomainLinearSystem> getSubdomainLinearSystem, Func<int, PsmSubdomainDofs> getSubdomainDofs)
		{
			this.environment = environment;
			this.subdomainTopology = subdomainTopology;
			this.getSubdomainLinearSystem = getSubdomainLinearSystem;
			this.getSubdomainDofs = getSubdomainDofs;
		}

		public IDictionary<int, DiagonalMatrix> SubdomainMatricesWb { get; } = new ConcurrentDictionary<int, DiagonalMatrix>();

		/// <summary>
		/// See eq (6.3) from Papagiannakis bachelor :
		/// Lpb^e = Db^e * Lb^e * inv( (Lb^e)^T * Db^e * Lb^e)
		/// </summary>
		public void CalcScalingMatrices(DistributedOverlappingIndexer indexer)
		{
			//Build Db^s from each subdomain's Kff
			Func<int, Vector> calcSubdomainDb = subdomainID =>
			{
				//ISubdomain subdomain = model.GetSubdomain(subdomainID);
				IMatrix Kff = getSubdomainLinearSystem(subdomainID).Matrix;

				//TODO: This should be a polymorphic method in the LinearAlgebra project. Interface IDiagonalizable 
				//		with methods: GetDiagonal() and GetSubdiagonal(int[] indices). Optimized versions for most storage
				//		formats are possible. E.g. for Symmetric CSR/CSC with ordered indices, the diagonal entry is the 
				//		last of each row/col. For general CSC/CSC with ordered indices, bisection can be used for to locate
				//		the diagonal entry of each row/col in log(nnzPerRow). In any case these should be hidden from DDM classes.
				//TODO: It would be better to extract the diagonal from Kbb directly.
				Vector Df = Kff.GetDiagonal();

				int[] boundaryDofs = getSubdomainDofs(subdomainID).DofsBoundaryToFree;
				Vector Db = Df.GetSubvector(boundaryDofs);

				return Db;
			};
			Dictionary<int, Vector> diagonalStiffnesses = environment.CalcNodeData(calcSubdomainDb);

			// Use distributed vectors to let each subdomain inform its neighbors about its stiffness at their common dofs
			var distributedVector = new DistributedOverlappingVector(indexer, diagonalStiffnesses);
			distributedVector.RegularizeOverlappingEntries();

			Action<int> storeRelativeStiffness = subdomainID =>
			{
				relativeStiffnesses[subdomainID] = distributedVector.LocalVectors[subdomainID].RawData;
				SubdomainMatricesWb[subdomainID] = DiagonalMatrix.CreateFromArray(relativeStiffnesses[subdomainID]);
			};
			environment.DoPerNode(storeRelativeStiffness);
		}

		public void ScaleBoundaryRhsVector(int subdomainID, Vector boundaryRhsVector)
		{
			double[] coefficients = relativeStiffnesses[subdomainID];
			Debug.Assert(boundaryRhsVector.Length == coefficients.Length);
			for (int i = 0; i < coefficients.Length; i++)
			{
				boundaryRhsVector[i] *= coefficients[i];
			}
		}

		private Dictionary<int, bool> MustComputeSubdomainData(IModifiedSubdomains modifiedSubdomains)
		{
			bool isFirstAnalysis = relativeStiffnesses.Count == 0;
			return environment.CalcNodeData(subdomainID =>
			{
				if (isFirstAnalysis || modifiedSubdomains.IsMatrixModified(subdomainID))
				{
					return true;
				}

				foreach (int neighborID in subdomainTopology.GetNeighborsOfSubdomain(subdomainID))
				{
					if (modifiedSubdomains.IsMatrixModified(neighborID))
					{
						return true;
					}
				}

				return false;
			});
		}
	}
}
