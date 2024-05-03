namespace MGroup.Solvers.DDM.LinearSystem
{
	public class NullModifiedSubdomains : IModifiedSubdomains
	{
		public bool IsConnectivityModified(int subdomainID) => true;

		public bool IsMatrixModified(int subdomainID) => true;

		public bool IsRhsModified(int subdomainID) => true;
	}
}
