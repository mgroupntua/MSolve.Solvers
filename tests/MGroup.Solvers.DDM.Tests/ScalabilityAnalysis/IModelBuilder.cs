namespace MGroup.Solvers.DDM.Tests.ScalabilityAnalysis
{
	using MGroup.Environments;
	using MGroup.MSolve.Discretization.Entities;
	using MGroup.Solvers.DDM.FetiDP.Dofs;

	public interface IModelBuilder
	{
		double[] DomainLengthPerAxis { get; set; }

		int[] NumElementsPerAxis { get; set; }

		int[] NumSubdomainsPerAxis { get; set; }

		int[] SubdomainSizePerElementSize { get; }

		(IModel model, ComputeNodeTopology nodeTopology) CreateMultiSubdomainModel();

		IModel CreateSingleSubdomainModel();

		ICornerDofSelection GetCornerDofs(IModel model);

		(List<int[]> numElements, int[] numSubdomains) GetParametricConfigConstNumSubdomains();

		(int[] numElements, List<int[]> numSubdomains) GetParametricConfigConstNumElements();

		(List<int[]> numElements, List<int[]> numSubdomains) GetParametricConfigConstSubdomainPerElementSize();
	}
}
