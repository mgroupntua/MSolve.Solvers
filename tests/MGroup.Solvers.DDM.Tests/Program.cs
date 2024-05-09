using MGroup.Solvers.DDM.Tests;

class Program
{
	static void Main(string[] args)
	{
		MpiTestSuite.RunTestsWith5Processes();
		//MpiScalabilityAnalysisRunner.RunScalabilityAnalysesWith4Processes();
	}
}
