using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Diagnostics;
using System.Text;
using MGroup.Environments;
using MGroup.MSolve.Discretization.Dofs;
using MGroup.Solvers.DDM.LinearSystem;

namespace MGroup.Solvers.DDM.FetiDP.Dofs
{
	public class GeneralModifiedCornerDofs : IModifiedCornerDofs
	{
		private readonly IComputeEnvironment environment;
		private readonly Func<int, FetiDPSubdomainDofs> getSubdomainDofs;

		private bool isFirstAnalysis;
		private readonly ConcurrentDictionary<int, IntDofTable> previousSubdomainCornerDofs 
			= new ConcurrentDictionary<int, IntDofTable>();
		private readonly ConcurrentDictionary<int, bool> subdomainsWithModifiedCornerDofs = new ConcurrentDictionary<int, bool>();

		public GeneralModifiedCornerDofs(IComputeEnvironment environment, Func<int, FetiDPSubdomainDofs> getSubdomainDofs)
		{
			this.environment = environment;
			this.getSubdomainDofs = getSubdomainDofs;
			isFirstAnalysis = true;
		}

		public bool AreGlobalCornerDofsModified { get; private set; }

		public bool AreSubdomainCornerDofsModified(int subdomainID) => subdomainsWithModifiedCornerDofs[subdomainID];

		public void Update(IModifiedSubdomains modifiedSubdomains)
		{
			if (isFirstAnalysis)
			{
				environment.DoPerNode(subdomainID =>
				{
					IntDofTable newCornerDofs = getSubdomainDofs(subdomainID).DofOrderingCorner;
					previousSubdomainCornerDofs[subdomainID] = newCornerDofs;
					subdomainsWithModifiedCornerDofs[subdomainID] = true;
					//LogMsg(subdomainID);
				});

				AreGlobalCornerDofsModified = true;
			}
			else
			{
				environment.DoPerNode(subdomainID =>
				{
					IntDofTable newCornerDofs = getSubdomainDofs(subdomainID).DofOrderingCorner;
					if (modifiedSubdomains.IsConnectivityModified(subdomainID))
					{
						IntDofTable previousCornerDofs = previousSubdomainCornerDofs[subdomainID];

						bool areModified = !previousCornerDofs.HasSameRowsColumns(newCornerDofs);
						previousSubdomainCornerDofs[subdomainID] = newCornerDofs;
						subdomainsWithModifiedCornerDofs[subdomainID] = areModified;
						//LogMsg(subdomainID);
					}
					else
					{
						previousSubdomainCornerDofs[subdomainID] = newCornerDofs;
						subdomainsWithModifiedCornerDofs[subdomainID] = false;
					}
				});

				AreGlobalCornerDofsModified = environment.AllReduceOr(subdomainsWithModifiedCornerDofs);

			}
			
			isFirstAnalysis = false;
			//LogGlobalMsg();
		}

		#region log
		//private void LogMsg(int subdomainID)
		//{
		//	if (subdomainsWithModifiedCornerDofs[subdomainID])
		//	{
		//		Console.WriteLine($"Corners dofs of subdomain {subdomainID} have been modified");
		//		Debug.WriteLine($"Corners dofs of subdomain {subdomainID} have been modified");
		//	}
		//}

		//private void LogGlobalMsg()
		//{
		//	if (AreGlobalCornerDofsModified)
		//	{
		//		Console.WriteLine("Global corners dofs have been modified");
		//		Debug.WriteLine("Global corners dofs have been modified");
		//	}
		//	else
		//	{
		//		Console.WriteLine("Global corners dofs are the same");
		//		Debug.WriteLine("Global corners dofs are the same");
		//	}
		//}
		#endregion
	}
}
