using MGroup.MSolve.Discretization.Entities;
using MGroup.Solvers.DDM.Commons;
using MGroup.Solvers.DDM.LinearSystem;

namespace MGroup.Solvers.DDM.FetiDP.Dofs
{
	public class FetiDPSubdomainDofs
	{
		private readonly ISubdomain subdomain;
		private readonly ISubdomainLinearSystem linearSystem;

		public FetiDPSubdomainDofs(ISubdomain subdomain, ISubdomainLinearSystem linearSystem)
		{
			this.subdomain = subdomain;
			this.linearSystem = linearSystem;
		}

		public IntDofTable DofOrderingCorner { get; private set; }

		public IntDofTable DofOrderingBoundaryRemainder { get; private set; }

		public int[] DofsBoundaryRemainderToRemainder { get; private set; }

		public int[] DofsCornerToFree { get; private set; }

		public int[] DofsInternalToRemainder { get; private set; }

		public int[] DofsRemainderToFree { get; private set; }

		public bool IsEmpty => DofOrderingCorner == null;

		public (int[] internalToFree, int[] boundaryRemainderToFree) RemapDofs()
		{
			var internalToFree = new int[DofsInternalToRemainder.Length];
			for (int i = 0; i < internalToFree.Length; ++i)
			{
				internalToFree[i] = DofsRemainderToFree[DofsInternalToRemainder[i]];
			}

			var boundaryRemainderToFree = new int[DofsBoundaryRemainderToRemainder.Length];
			for (int i = 0; i < boundaryRemainderToFree.Length; ++i)
			{
				boundaryRemainderToFree[i] = DofsRemainderToFree[DofsBoundaryRemainderToRemainder[i]];
			}

			return (internalToFree, boundaryRemainderToFree);
		}

		public void ReorderInternalDofs(DofPermutation permutation)
		{
			if (permutation.IsBetter)
			{
				DofsInternalToRemainder = permutation.ReorderKeysOfDofIndicesMap(DofsInternalToRemainder);
			}
		}

		public void ReorderRemainderDofs(DofPermutation permutation)
		{
			if (permutation.IsBetter)
			{
				DofsRemainderToFree = permutation.ReorderKeysOfDofIndicesMap(DofsRemainderToFree);

				if (DofsInternalToRemainder != null)
				{
					DofsInternalToRemainder = permutation.ReorderValuesOfDofIndicesMap(DofsInternalToRemainder);
				}
				if (DofsBoundaryRemainderToRemainder != null)
				{
					DofsBoundaryRemainderToRemainder = permutation.ReorderValuesOfDofIndicesMap(DofsBoundaryRemainderToRemainder);
				}
			}
		}

		public void SeparateAllFreeDofs(ICornerDofSelection cornerDofSelection)
		{
			var cornerDofOrdering = new IntDofTable();
			var cornerToFree = new List<int>();
			var remainderToFree = new List<int>();
			var boundaryRemainderDofOrdering = new IntDofTable();
			var boundaryRemainderToRemainder = new List<int>();
			var internalToRemainder = new List<int>();
			IntDofTable freeDofs = linearSystem.DofOrdering.FreeDofs;
			IEnumerable<int> nodes = freeDofs.GetRows(); //TODO: Optimize access: Directly get INode, Dictionary<IDof, int>
			foreach (int nodeID in nodes)
			{
				int nodeMultiplicity = subdomain.GetMultiplicityOfNode(nodeID);
				IReadOnlyDictionary<int, int> freeDofsOfNode = freeDofs.GetDataOfRow(nodeID);
				foreach (var dofIdxPair in freeDofsOfNode)
				{
					int dofID = dofIdxPair.Key;
					int freeDofIdx = dofIdxPair.Value;
					if (cornerDofSelection.IsCornerDof(subdomain.ID, nodeID, dofID))
					{
						cornerDofOrdering[nodeID, dofID] = cornerToFree.Count;
						cornerToFree.Add(freeDofIdx);
					}
					else
					{
						int remainderDofIdx = remainderToFree.Count;
						remainderToFree.Add(freeDofIdx);

						if (nodeMultiplicity > 1)
						{
							// Boundary remainder dof
							boundaryRemainderDofOrdering[nodeID, dofID] = boundaryRemainderToRemainder.Count;
							boundaryRemainderToRemainder.Add(remainderDofIdx);
						}
						else
						{
							// Internal dof
							internalToRemainder.Add(remainderDofIdx);
						}
					}
				}
			}

			this.DofOrderingCorner = cornerDofOrdering;
			this.DofsCornerToFree = cornerToFree.ToArray();
			this.DofsRemainderToFree = remainderToFree.ToArray();
			this.DofOrderingBoundaryRemainder = boundaryRemainderDofOrdering;
			this.DofsBoundaryRemainderToRemainder = boundaryRemainderToRemainder.ToArray();
			this.DofsInternalToRemainder = internalToRemainder.ToArray();
		}

		public void SeparateFreeDofsIntoCornerAndRemainder(ICornerDofSelection cornerDofSelection)
		{
			var cornerDofOrdering = new IntDofTable();
			var cornerToFree = new List<int>();
			var remainderToFree = new HashSet<int>();
			int numCornerDofs = 0;
			IntDofTable freeDofs = linearSystem.DofOrdering.FreeDofs;
			IEnumerable<int> nodes = freeDofs.GetRows(); //TODO: Optimize access: Directly get INode, Dictionary<IDof, int>
			foreach (int node in nodes)
			{
				IReadOnlyDictionary<int, int> dofsOfNode = freeDofs.GetDataOfRow(node);
				foreach (var dofIdxPair in dofsOfNode)
				{
					int dof = dofIdxPair.Key;
					if (cornerDofSelection.IsCornerDof(subdomain.ID, node, dof))
					{
						cornerDofOrdering[node, dof] = numCornerDofs++;
						cornerToFree.Add(dofIdxPair.Value);
					}
					else
					{
						remainderToFree.Add(dofIdxPair.Value);
					}
				}
			}

			this.DofOrderingCorner = cornerDofOrdering;
			this.DofsCornerToFree = cornerToFree.ToArray();
			this.DofsRemainderToFree = remainderToFree.ToArray();
		}
	}
}
