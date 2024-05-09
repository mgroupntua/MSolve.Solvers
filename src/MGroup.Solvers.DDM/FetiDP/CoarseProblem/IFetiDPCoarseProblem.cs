namespace MGroup.Solvers.DDM.FetiDP.CoarseProblem
{
	using MGroup.LinearAlgebra.Vectors;
	using MGroup.Solvers.DDM.FetiDP.Dofs;
	using MGroup.Solvers.DDM.Output;

	public interface IFetiDPCoarseProblem
	{
		void FindCoarseProblemDofs(DdmLogger logger, IModifiedCornerDofs modifiedCornerDofs);

		void PrepareMatricesForSolution();

		//TODOMPI: Perhaps I should have dedicated classes for these distributed vectors. DistributedOverlappingVector was 
		//		cumbersome, because it required an indexer. Having an indexer for coarse-problem vectors makes sense only in 
		//		distributed implementations. Right now, in distributed implementations, these dictionaries contain vectors
		//		that can be local to different memory spaces, which might be confusing.
		void SolveCoarseProblem(IDictionary<int, Vector> coarseProblemRhs, IDictionary<int, Vector> coarseProblemSolution);
	}
}
