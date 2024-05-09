
namespace MGroup.Solvers.DDM.FetiDP.Dofs
{
	public class UserDefinedCornerDofSelection : ICornerDofSelection
	{
		private readonly HashSet<int> cornerNodes = new HashSet<int>();

		public UserDefinedCornerDofSelection()
		{
		}

		public int[] CornerNodeIDs => cornerNodes.ToArray();

		public void AddCornerNode(int nodeID) => cornerNodes.Add(nodeID);

		public bool IsCornerDof(int subdomainID, int nodeID, int dofID) => cornerNodes.Contains(nodeID);
	}
}
