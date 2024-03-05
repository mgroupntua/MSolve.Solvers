using System;
using System.Collections.Generic;
using System.Text;
using MGroup.MSolve.Discretization;

namespace MGroup.Solvers.DDM.LagrangeMultipliers
{
	public interface ICrossPointStrategy
	{
		List<(int subdomainPlus, int subdomainMinus)> ListSubdomainCombinations(IEnumerable<int> subdomainIDs);
	}
}
