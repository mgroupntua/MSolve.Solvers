using MGroup.Solvers.DDM.LinearSystem;

namespace MGroup.Solvers.DDM.FetiDP.Dofs
{
	public interface IModifiedCornerDofs
	{
		bool AreGlobalCornerDofsModified { get; }

		bool AreSubdomainCornerDofsModified(int subdomainID);

		void Update(IModifiedSubdomains modifiedSubdomains);
	}
}
