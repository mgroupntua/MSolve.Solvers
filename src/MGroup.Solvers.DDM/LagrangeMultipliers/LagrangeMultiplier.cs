using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Text;

namespace MGroup.Solvers.DDM.LagrangeMultipliers
{
	public class LagrangeMultiplier : IComparable<LagrangeMultiplier>
	{
		public LagrangeMultiplier(int nodeID, int dofID, int subdomainPlus, int subdomainMinus, int localIdx)
		{
			NodeID = nodeID;
			DofID = dofID;
			SubdomainMinus = subdomainMinus;
			SubdomainPlus = subdomainPlus;
			this.LocalIdx = localIdx;
		}

		public int DofID { get; }

		public int LocalIdx { get; }

		public int NodeID { get; }

		public int SubdomainMinus { get; }

		public int SubdomainPlus { get; }

		public int CompareTo(LagrangeMultiplier other)
		{
			if (this.NodeID - other.NodeID != 0)
			{
				return this.NodeID - other.NodeID;
			}
			else if (this.DofID - other.DofID != 0)
			{
				return this.DofID - other.DofID;
			}
			else if ((this.SubdomainPlus != other.SubdomainPlus) || (this.SubdomainMinus != other.SubdomainMinus))
			{
				throw new NotImplementedException("We cannot compare lagranges between the same dofs, " +
					"but different subdomains at crosspoints yet.");
			}
			else
			{
				return 0;
			}
		}
	}
}
