namespace MGroup.Solvers.DDM.FetiDP.Dofs
{
	public interface ICornerDofSelection
	{
		bool IsCornerDof(int subdomainID, int nodeID, int dofID);
	}
}
