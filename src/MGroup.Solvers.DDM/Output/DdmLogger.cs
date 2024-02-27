using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Text;
using MGroup.Environments;
using MGroup.Solvers.Logging;

namespace MGroup.Solvers.DDM.Output
{
	public class DdmLogger
	{
		private readonly IComputeEnvironment environment;
		private readonly string solverName;
		private readonly int numSubdomains;
		private int analysisIteration = -1;

		private SortedDictionary<int, (int numIterations, double residualNormRatio)> convergenceData
			= new SortedDictionary<int, (int numIterations, double residualNormRatio)>();
		private List<Dictionary<int, int>> problemSizeData = new List<Dictionary<int, int>>();
		private List<Dictionary<int, int>> subdomainProblemSizeData = new List<Dictionary<int, int>>();
		private List<(int local, int remote)> totalTransfersData = new List<(int local, int remote)>();

		public DdmLogger(IComputeEnvironment environment, string solverName, int numSubdomains)
		{
			this.environment = environment;
			this.solverName = solverName;
			this.numSubdomains = numSubdomains;
		}

		public string Header { get; set; }

		public void IncrementAnalysisIteration()
		{
			++analysisIteration;
			problemSizeData.Add(new Dictionary<int, int>());
			subdomainProblemSizeData.Add(new Dictionary<int, int>());
			totalTransfersData.Add((-1, -1));
		}

		/// <summary>
		/// Log the number of unique dofs for a specified level: Lvl 0 = free dofs, Lvl 1 = interface problem, 
		/// Lvl 2 = coarse problem, Lvl 3 = interface problem of coarse problem (multilevel-DDM), 
		/// Lvl 4 = coarse problem of coarse problem (multilevel-DDM), ...
		/// </summary>
		/// <param name="problemLevel"></param>
		/// <param name="size"></param>
		public void LogProblemSize(int problemLevel, int size)
		{
			problemSizeData[analysisIteration][problemLevel] = size;
		}

		public void LogSolverConvergenceData(int numIterations, double residualNormRatio)
		{
			convergenceData[analysisIteration] = (numIterations, residualNormRatio);
		}

		public void LogSubdomainProblemSize(int subdomainID, int subdomainSize)
		{
			subdomainProblemSizeData[analysisIteration][subdomainID] = subdomainSize;
		}

		public void LogTransfers(int totalLocalTransfers, int totalRemoteTransfers)
		{
			totalTransfersData[analysisIteration] = (totalLocalTransfers, totalRemoteTransfers);
		}

		public void WriteToConsole()
		{
			environment.DoGlobalOperation(() =>
			{
				// Point the stream to standard output
				var writer = new StreamWriter(Console.OpenStandardOutput());
				writer.AutoFlush = true;
				Console.SetOut(writer);

				// Call the abstract method
				WriteToStream(writer);

				// Recover the standard output stream
				var standardOutput = new StreamWriter(Console.OpenStandardOutput());
				standardOutput.AutoFlush = true;
				Console.SetOut(standardOutput);
			});
		}

		public void WriteToDebug()
		{
			throw new NotImplementedException();
		}

		public void WriteToFile(string path, bool append = true)
		{
			environment.DoGlobalOperation(() =>
			{
				using (var writer = new StreamWriter(path, append))
				{
#if DEBUG
					writer.AutoFlush = true; // To look at intermediate output at certain breakpoints
#endif
					WriteToStream(writer);
				}
			});
		}

		public void WriteToStream(StreamWriter writer)
		{
			writer.WriteLine("************************************************************************************************");
			WriteHeader(writer);
			for (int t = 0; t < convergenceData.Count; ++t)
			{
				(int numIterations, double residualNorm) = convergenceData[t];
				var msgConvergence = new StringBuilder();
				msgConvergence.Append($"Analysis iteration {t}:");
				msgConvergence.Append($" Solver iterations = {numIterations}.");
				msgConvergence.Append($" Residual norm ratio = {residualNorm}.");
				writer.WriteLine(msgConvergence);

				var msgSizes = new StringBuilder();
				msgSizes.Append($"Analysis iteration {t}:");
				if (problemSizeData[t].TryGetValue(0, out int size0))
				{
					msgSizes.Append($" Global problem size = {size0}.");
				}
				if (problemSizeData[t].TryGetValue(1, out int size1))
				{
					msgSizes.Append($" Interface problem size = {size1}.");
				}
				if (problemSizeData[t].TryGetValue(2, out int size2))
				{
					msgSizes.Append($" Coarse problem size = {size2}.");
				}
				msgSizes.AppendLine();

				try
				{
					(int subSizeMin, int subSizeMax, _, double subSizeAvg, double subSizeStdev) =
						LoggingUtilities.CalcStatistics(subdomainProblemSizeData[t].Values);
					msgSizes.AppendLine(
						$"\t Subdomain problem size: min={subSizeMin}, max={subSizeMax}, avg={subSizeAvg}, stdev={subSizeStdev}.");

					(int localTransfers, int remoteTransfers) = totalTransfersData[t];
					msgSizes.AppendLine(
						$"\t Transfers per global vector operation: local={localTransfers}, remote={remoteTransfers}");

					writer.Write(msgSizes);
				}
				catch (Exception)
				{
				}
			}

			(int iterMin, int iterMax, double iterAvg, double iterSum, double iterStdev) 
				= LoggingUtilities.CalcStatistics(convergenceData.Select(d => d.Value.numIterations));
			writer.WriteLine($"Solver iteration statistics: min = {iterMin}, max = {iterMax}, sum = {iterSum}, " +
				$"average = {iterAvg}, stdev = {iterStdev}");
			try
			{
				WriteAggregates(writer);
			}
			catch (Exception)
			{
			}
			
			writer.WriteLine("************************************************************************************************");
		}

		private void WriteAggregates(StreamWriter writer)
		{
			(int globalDofsMin, int globalDofsMax, _, double globalDofsAvg, double globalDofsStdev)
					= LoggingUtilities.CalcStatistics(problemSizeData.Select(d => d[0]));
			writer.WriteLine($"Global problem size statistics: min = {globalDofsMin}, max = {globalDofsMax}, " +
				$"average = {globalDofsAvg}, stdev = {globalDofsStdev}");

			(int interfaceDofsMin, int interfaceDofsMax, _, double interfaceDofsAvg, double interfaceDofsStdev)
				= LoggingUtilities.CalcStatistics(problemSizeData.Select(d => d[1]));
			writer.WriteLine($"Interface problem size statistics: min = {interfaceDofsMin}, max = {interfaceDofsMax}, " +
				$"average = {interfaceDofsAvg}, stdev = {interfaceDofsStdev}");

			(int coarseDofsMin, int coarseDofsMax, _, double coarseDofsAvg, double coarseDofsStdev)
				= LoggingUtilities.CalcStatistics(problemSizeData.Select(d => d[2]));
			writer.WriteLine($"Coarse problem size statistics: min = {coarseDofsMin}, max = {coarseDofsMax}, " +
				$"average = {coarseDofsAvg}, stdev = {coarseDofsStdev}");


			try
			{
				var allSubdomainSizes = new List<int>(1 + numSubdomains * subdomainProblemSizeData.Count);
				foreach (Dictionary<int, int> sizePerSubdomain in subdomainProblemSizeData)
				{
					allSubdomainSizes.AddRange(sizePerSubdomain.Values);
				}

				(int subdomainDofsMin, int subdomainDofsMax, _, double subdomainDofsAvg, double subdomainDofsStdev)
					= LoggingUtilities.CalcStatistics(allSubdomainSizes);
				writer.WriteLine($"Subdomain problem size statistics: min = {subdomainDofsMin}, max = {subdomainDofsMax}, " +
					$"average = {subdomainDofsAvg}, stdev = {subdomainDofsStdev}");

				(int localTransfersMin, int localTransfersMax, _, double localTransfersAvg, double localTransfersStdev)
					= LoggingUtilities.CalcStatistics(totalTransfersData.Select(t => t.local));
				writer.WriteLine($"Local transfers per global vector operation statistics: min = {localTransfersMin}, " +
					$"max = {localTransfersMax}, average = {localTransfersAvg}, stdev = {localTransfersStdev}");

				(int remoteTransfersMin, int remoteTransfersMax, _, double remoteTransfersAvg, double remoteTransfersStdev)
					= LoggingUtilities.CalcStatistics(totalTransfersData.Select(t => t.remote));
				writer.WriteLine($"Remote transfers per global vector operation statistics: min = {remoteTransfersMin}, " +
					$"max = {remoteTransfersMax}, average = {remoteTransfersAvg}, stdev = {remoteTransfersStdev}");
			}
			catch (Exception)
			{
			}
		}

		private void WriteHeader(StreamWriter writer)
		{
			if (Header != null)
			{
				writer.Write(Header);
			}
			else
			{
				writer.WriteLine(DateTime.Now);
				writer.WriteLine($"Solver: {solverName}");
				writer.WriteLine($"Num subdomains: {numSubdomains}");
			}
		}
	}
}
