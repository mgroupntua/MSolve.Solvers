namespace MGroup.Solvers.DDM.FetiDP.CoarseProblem
{
	using MGroup.Environments;
	using MGroup.LinearAlgebra.Matrices;
	using MGroup.LinearAlgebra.Vectors;
	using MGroup.Solvers.DDM.Commons;
	using MGroup.Solvers.DDM.FetiDP.Dofs;
	using MGroup.Solvers.DDM.FetiDP.StiffnessMatrices;
	using MGroup.Solvers.DDM.Output;

	public class FetiDPCoarseProblemGlobal : IFetiDPCoarseProblem
	{
		private readonly IComputeEnvironment environment;
		private readonly FetiDPCoarseProblemGlobalDofs coarseProblemDofs;
		private readonly IFetiDPCoarseProblemGlobalMatrix coarseProblemMatrix;
		private readonly FetiDPCoarseProblemGlobalSolver coarseProblemSolver;
		private readonly Func<int, FetiDPSubdomainDofs> getSubdomainDofs;
		private readonly Func<int, IFetiDPSubdomainMatrixManager> getSubdomainMatrices;

		public FetiDPCoarseProblemGlobal(IComputeEnvironment environment, IFetiDPCoarseProblemGlobalMatrix coarseProblemMatrix,
			Func<int, FetiDPSubdomainDofs> getSubdomainDofs, Func<int, IFetiDPSubdomainMatrixManager> getSubdomainMatrices)
		{
			this.environment = environment;
			this.coarseProblemMatrix = coarseProblemMatrix;
			this.getSubdomainDofs = getSubdomainDofs;
			this.getSubdomainMatrices = getSubdomainMatrices;
			this.coarseProblemDofs = new FetiDPCoarseProblemGlobalDofs();
			this.coarseProblemSolver = new FetiDPCoarseProblemGlobalSolver(coarseProblemDofs, coarseProblemMatrix);
		}

		public void FindCoarseProblemDofs(DdmLogger logger, IModifiedCornerDofs modifiedCornerDofs)
		{
			if (!modifiedCornerDofs.AreGlobalCornerDofsModified)
			{
				#region log
				//Console.WriteLine("Coarse problem dofs are the same as last analysis");
				//Debug.WriteLine("Coarse problem dofs are the same as last analysis");
				#endregion

				if (logger != null)
				{
					logger.LogProblemSize(2, coarseProblemDofs.NumGlobalCornerDofs);
				}
				return;
			}

			Dictionary<int, IntDofTable> subdomainCornerDofs =
				environment.CalcNodeDataAndTransferToGlobalMemory(s => getSubdomainDofs(s).DofOrderingCorner);
			//Dictionary<int, IntDofTable> subdomainCornerDofs = RegatherSubdomainCornerDofs(modifiedCornerDofs);

			environment.DoGlobalOperation(() =>
			{
				IEnumerable<int> subdomainIDs = subdomainCornerDofs.Keys;
				coarseProblemDofs.FindGlobalCornerDofs(subdomainCornerDofs);
				coarseProblemDofs.SubdomainToGlobalCornerDofs = environment.DoPerItemInGlobalMemory(
					subdomainIDs, s => coarseProblemDofs.CalcSubdomainGlobalCornerDofMap(s));

				DofPermutation permutation = coarseProblemMatrix.ReorderGlobalCornerDofs(
					coarseProblemDofs.NumGlobalCornerDofs, coarseProblemDofs.SubdomainToGlobalCornerDofs);
				if (permutation.IsBetter)
				{
					// The global dof ordering and the subdomain-global maps need to be updated.
					coarseProblemDofs.GlobalDofOrderingCorner.Reorder(
						permutation.PermutationArray, permutation.PermutationIsOldToNew);
					coarseProblemDofs.SubdomainToGlobalCornerDofs = environment.DoPerItemInGlobalMemory(
						subdomainIDs, s => coarseProblemDofs.CalcSubdomainGlobalCornerDofMap(s));
				}

				if (logger != null)
				{
					logger.LogProblemSize(2, coarseProblemDofs.NumGlobalCornerDofs);
				}
			});
		}

		public void PrepareMatricesForSolution()
		{
			// This is done by the solver itself
			//environment.DoPerNode(subdomainID => getSubdomainMatrices(subdomainID).CalcSchurComplementOfRemainderDofs());

			Dictionary<int, IMatrix> subdomainMatricesScc =
				environment.CalcNodeDataAndTransferToGlobalMemory(s => getSubdomainMatrices(s).SchurComplementOfRemainderDofs);

			environment.DoGlobalOperation(() =>
			{
				coarseProblemMatrix.Clear();
				coarseProblemMatrix.InvertGlobalScc(
					coarseProblemDofs.NumGlobalCornerDofs, coarseProblemDofs.SubdomainToGlobalCornerDofs, subdomainMatricesScc);
			});
		}

		public void SolveCoarseProblem(
			IDictionary<int, Vector> coarseProblemRhs, IDictionary<int, Vector> coarseProblemSolution)
		{
			Dictionary<int, Vector> subdomainRhsGlobal = 
				environment.CalcNodeDataAndTransferToGlobalMemory(s => coarseProblemRhs[s]);

			Vector globalSolution = null;
			environment.DoGlobalOperation(() =>
				globalSolution = coarseProblemSolver.SolveCoarseProblem(subdomainRhsGlobal)
			);

			//TODOMPI: write directly into vectors of coarseProblemSolution
			Dictionary<int, Vector> subdomainSolutionsLocal = environment.CalcNodeDataAndTransferToLocalMemory(
				s => coarseProblemSolver.ExtractCoarseProblemSolutionForSubdomain(s, globalSolution));

			environment.DoPerNode(s => coarseProblemSolution[s].CopyFrom(subdomainSolutionsLocal[s]));
		}

		/// <summary>
		/// Avoids transferring corner dofs from subdomains, if they are already available.
		/// </summary>
		/// <param name="modifiedCornerDofs"></param>
		private Dictionary<int, IntDofTable> RegatherSubdomainCornerDofs(IModifiedCornerDofs modifiedCornerDofs)
		{
			//TODO: This will probably not improve performance. Perhaps for MPI environments.
			//TODO: Changing the order of entries in the returned dictionary, seems to alter the solution. Find out why.

			if ((modifiedCornerDofs is NullModifiedCornerDofs) || (coarseProblemDofs.SubdomainDofOrderingsCorner == null))
			{
				return environment.CalcNodeDataAndTransferToGlobalMemory(s => getSubdomainDofs(s).DofOrderingCorner);
			}
			else
			{
				Dictionary<int, IntDofTable> transferedDofs = environment.CalcNodeDataAndTransferToGlobalMemoryPartial(
				   s => getSubdomainDofs(s).DofOrderingCorner, s => modifiedCornerDofs.AreSubdomainCornerDofsModified(s));

				#region log
				//foreach (int s in transferedDofs.Keys)
				//{
				//	Console.WriteLine($"Transfered corner dofs of subdomain {s} to global memory");
				//	Debug.WriteLine($"Transfered corner dofs of subdomain {s} to global memory");
				//}
				#endregion

				var result = new Dictionary<int, IntDofTable>(coarseProblemDofs.SubdomainDofOrderingsCorner.Count);
				foreach (KeyValuePair<int, IntDofTable> pair in coarseProblemDofs.SubdomainDofOrderingsCorner)
				{
					int subdomainID = pair.Key;
					bool isModified = transferedDofs.TryGetValue(subdomainID, out IntDofTable subdomainCornerDofs);
					result[subdomainID] = isModified ? subdomainCornerDofs : pair.Value;
				}

				return result;
			}
		}

		public class Factory : IFetiDPCoarseProblemFactory
		{
			private readonly IFetiDPCoarseProblemGlobalMatrix coarseProblemMatrix;

			public Factory(IFetiDPCoarseProblemGlobalMatrix coarseProblemMatrix)
			{
				this.coarseProblemMatrix = coarseProblemMatrix;
			}

			public IFetiDPCoarseProblem CreateCoarseProblem(
				IComputeEnvironment environment, ISubdomainTopology subdomainTopology,
				Func<int, FetiDPSubdomainDofs> getSubdomainDofs, Func<int, IFetiDPSubdomainMatrixManager> getSubdomainMatrices)
			{
				return new FetiDPCoarseProblemGlobal(environment, coarseProblemMatrix, getSubdomainDofs, getSubdomainMatrices);
			}
		}
	}
}
