using System;
using System.Collections.Generic;
using System.Text;

namespace MGroup.Solvers.Logging
{
	public static class LoggingUtilities
	{
		public static (int min, int max, int sum, double avg, double stdev) CalcStatistics(IEnumerable<int> sample)
		{
			int min = int.MaxValue;
			int max = int.MinValue;
			int sum = 0;
			int count = 0;
			foreach (int val in sample)
			{
				if (val < min)
				{
					min = val;
				}
				if (val > max)
				{
					max = val;
				}
				sum += val;
				++count;
			}
			double avg = ((double)sum) / count;

			double sumDev = 0;
			foreach (int val in sample)
			{
				sumDev += (val - avg) * (val - avg);
			}
			double stdev = Math.Sqrt(sumDev / (count - 1));

			return (min, max, sum, avg, stdev);
		}

		public static (double min, double max, double sum, double avg, double stdev) CalcStatistics(IEnumerable<double> sample)
		{
			double min = double.MaxValue;
			double max = double.MinValue;
			double sum = 0.0;
			int count = 0;
			foreach (double val in sample)
			{
				if (val < min)
				{
					min = val;
				}
				if (val > max)
				{
					max = val;
				}
				sum += val;
				++count;
			}
			double avg = sum / count;

			double sumDev = 0;
			foreach (double val in sample)
			{
				sumDev += (val - avg) * (val - avg);
			}
			double stdev = Math.Sqrt(sumDev / (count - 1));

			return (min, max, sum, avg, stdev);
		}
	}
}
