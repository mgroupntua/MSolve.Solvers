using System.Diagnostics;
using MGroup.Solvers.DDM.FetiDP.Dofs;
using MGroup.Solvers.DDM.Mappings;
using MGroup.Solvers.DDM.PSM.Dofs;

namespace MGroup.Solvers.DDM.PFetiDP.Dofs
{
	public class PFetiDPSubdomainDofs
	{
		private readonly PsmSubdomainDofs psmDofs;
		private readonly FetiDPSubdomainDofs fetiDPDofs;

		public PFetiDPSubdomainDofs(PsmSubdomainDofs psmDofs, FetiDPSubdomainDofs fetiDPDofs)
		{
			this.psmDofs = psmDofs;
			this.fetiDPDofs = fetiDPDofs;
		}

		public bool IsEmpty => MatrixNcb == null;


		/// <summary>
		/// Boolean mapping matrix where rows = corner dofs of subdomain, columns = boundary dofs of subdomain.
		/// </summary>
		public IMappingMatrix MatrixNcb { get; private set; }

		/// <summary>
		/// Boolean mapping matrix where rows = remainder dofs of subdomain, columns = boundary (not boundary remainder) dofs of 
		/// subdomain.
		/// </summary>
		public IMappingMatrix MatrixNrb { get; private set; }

		public void MapPsmFetiDPDofs()
		{
			// Free to boundary dofs
			int[] boundaryToFree = psmDofs.DofsBoundaryToFree;
			var freeToBoundary = new Dictionary<int, int>();
			for (int i = 0; i < boundaryToFree.Length; i++)
			{
				freeToBoundary[boundaryToFree[i]] = i;
			}

			// Corner to boundary dofs
			int[] cornerToFree = fetiDPDofs.DofsCornerToFree;
			var cornerToBoundary = new int[cornerToFree.Length];
			for (int c = 0; c < cornerToFree.Length; c++)
			{
				int f = cornerToFree[c];
				bool exists = freeToBoundary.TryGetValue(f, out int b); // all corner dofs are also boundary.
				Debug.Assert(exists, "Found corner dof that is not boundary. This should not have happened");
				cornerToBoundary[c] = b;
			}
			this.MatrixNcb = new BooleanMatrixRowsToColumns(cornerToFree.Length, boundaryToFree.Length, cornerToBoundary);

			// Remainder to boundary dofs
			int[] remainderToFree = fetiDPDofs.DofsRemainderToFree;
			var remainderToBoundary = new Dictionary<int, int>();
			for (int r = 0; r < remainderToFree.Length; r++)
			{
				int f = remainderToFree[r];
				bool exists = freeToBoundary.TryGetValue(f, out int b);
				if (exists) // some remainder dofs are internal, thus they cannot be boundary too.
				{
					remainderToBoundary[r] = b;
				}
			}
			this.MatrixNrb = new MappingMatrixN(remainderToFree.Length, boundaryToFree.Length, remainderToBoundary);
		}
	}
}
