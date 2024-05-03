namespace MGroup.Solvers.DDM.FetiDP.CoarseProblem
{
	using MGroup.Environments;
	using MGroup.Solvers.DDM.FetiDP.Dofs;
	using MGroup.Solvers.DDM.FetiDP.StiffnessMatrices;

	public interface IFetiDPCoarseProblemFactory
	{
		IFetiDPCoarseProblem CreateCoarseProblem(
			IComputeEnvironment environment, ISubdomainTopology subdomainTopology,
			Func<int, FetiDPSubdomainDofs> getSubdomainDofs, Func<int, IFetiDPSubdomainMatrixManager> getSubdomainMatrices);
	}
}
