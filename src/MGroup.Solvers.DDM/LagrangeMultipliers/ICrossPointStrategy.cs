namespace MGroup.Solvers.DDM.LagrangeMultipliers
{
	using System.Collections.Generic;

	public interface ICrossPointStrategy
	{
		List<(int subdomainPlus, int subdomainMinus)> ListSubdomainCombinations(IEnumerable<int> subdomainIDs);
	}
}
