using MGroup.LinearAlgebra.Vectors;

namespace MGroup.Solvers.DDM.FetiDP.CoarseProblem
{
	public class FetiDPCoarseProblemGlobalSolver
	{
		private readonly IFetiDPCoarseProblemGlobalMatrix coarseProblemMatrix;
		private readonly FetiDPCoarseProblemGlobalDofs coarseProblemDofs;

		public FetiDPCoarseProblemGlobalSolver(FetiDPCoarseProblemGlobalDofs coarseProblemDofs,
			IFetiDPCoarseProblemGlobalMatrix coarseProblemMatrix)
		{
			this.coarseProblemDofs = coarseProblemDofs;
			this.coarseProblemMatrix = coarseProblemMatrix;
		}

		public Vector SolveCoarseProblem(IDictionary<int, Vector> coarseProblemRhs)
		{
			// Map reduce subdomain vectors to global
			var globalRhs = Vector.CreateZero(coarseProblemDofs.NumGlobalCornerDofs);
			foreach (int s in coarseProblemRhs.Keys)
			{
				Vector subdomainRhs = coarseProblemRhs[s];

				int[] subdomainToGlobalDofs = coarseProblemDofs.SubdomainToGlobalCornerDofs[s];
				globalRhs.AddIntoThisNonContiguouslyFrom(subdomainToGlobalDofs, subdomainRhs);

				//TODOMPI: delete this and the Lc matrices. They are used below, but this is not faster than the code above.
				//Mappings.IMappingMatrix Lc = coarseProblemDofs.SubdomainMatricesLc[s];
				//globalRhs.AddIntoThis(Lc.Multiply(subdomainRhs, true));
			}

			// Solve global problem
			var globalSolution = Vector.CreateZero(coarseProblemDofs.NumGlobalCornerDofs);
			coarseProblemMatrix.MultiplyInverseScc(globalRhs, globalSolution);

			return globalSolution;
		}

		public Vector ExtractCoarseProblemSolutionForSubdomain(int subdomainID, Vector globalSolution)
		{
			int[] subdomainToGlobalDofs = coarseProblemDofs.SubdomainToGlobalCornerDofs[subdomainID];
			Vector subdomainSolution = globalSolution.GetSubvector(subdomainToGlobalDofs);
			return subdomainSolution; //TODO: write into the existing vector directly during GetSubvector()

			//TODO: delete this and the Lc matrices. They are used below, but this is not faster than the code above.
			//Mappings.IMappingMatrix Lc = coarseProblemDofs.SubdomainMatricesLc[subdomainID];
			//Vector subvector = Lc.Multiply (globalSolution, false);
			//subdomainSolution.CopyFrom(subvector); //TODO: write into the existing vector directly during GetSubvector()
		}
	}
}
