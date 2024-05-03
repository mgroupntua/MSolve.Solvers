//TODO: Needs a much better name.
namespace MGroup.Solvers.DDM.LinearSystem
{
	public interface IModifiedSubdomains
	{
		bool IsConnectivityModified(int subdomainID);

		bool IsMatrixModified(int subdomainID);

		bool IsRhsModified(int subdomainID);
	}
}
