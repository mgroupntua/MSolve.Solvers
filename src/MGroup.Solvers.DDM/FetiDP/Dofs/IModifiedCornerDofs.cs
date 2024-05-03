namespace MGroup.Solvers.DDM.FetiDP.Dofs
{
	using MGroup.Solvers.DDM.LinearSystem;

	public interface IModifiedCornerDofs
	{
		bool AreGlobalCornerDofsModified { get; }

		bool AreSubdomainCornerDofsModified(int subdomainID);

		void Update(IModifiedSubdomains modifiedSubdomains);
	}
}
